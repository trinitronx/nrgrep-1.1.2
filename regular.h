
/*

Nrgrep -- a fast and flexible pattern matching tool.
Copyright (C) 2000,2001 Gonzalo Navarro

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of 
Chile. Blanco Encalada 2120, Santiago, Chile. gnavarro@dcc.uchile.cl

*/

#ifndef REGULARINCLUDED
#define REGULARINCLUDED

	/* Search a regular expression in a text */

#include "basics.h"	/* basics */
#include "options.h"	/* buffer management */
#include "parser.h"	/* parser */

typedef struct
   { mask **dtrans;     /* deterministic reduced automaton */
     mask dinitial;     /* initial states in reduced automaton */
     mask dfinal;       /* final states in reduced automaton */
     mask dtransin[256];/* precomputed for skip loop */
     mask B[256];	/* B table for reduced automaton */
     int dm;            /* number of states in reduced automaton */
     int slices;        /* number of deterministic slices */
     int width;         /* number of bits per slice */
     int wlen;          /* window length for bwd scan */
   } regularScanData;

typedef struct 

   { scanProc scanText; /* procedure for fast scanning */
     freeScanProc scanFree; /* procedure for freeing scanText */
     void *scanData;    /* data for fast scanning */
     int m,dm;             /* number of states, and of filtering one */
     int type;		/* FWD or BWD */
     mask match;	/* reduced states that matched */

     int width;         /* shape of verif autom */
     Mask initial,final;/* of fwd verif autom */
     Mask current;      /* state of verif */
     Mask **dftransV;   /* deterministic version of fwd verif automaton */
     Mask **dbtransV;   /* deterministic version of bwd verif automaton */
     Mask B[256];       /* B masks for verification */
     int *revMap;	/* map local -> global states */
     Mask V1;		 /* extra space for verification */
   } regularData;

	/* Preprocesses pat and creates a regularData structure
	   for searching */

regularData *regularPreproc (byte *pat, Tree *tree, Tree **pos);

	/* Frees P */

void regularFree (regularData *P);

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

bool regularSearch (byte **beg, byte **end, regularData *P);

        /* compute number of bits necessary to hold the regexp */

void regularLength (Tree *e, int *L);

        /* remaps the states so that only the active ones exist. it does
           not check that the reduction is consistent. width is the width
           that will be used for the deterministic tables (<= W) and the
           mapping guarantees that no slice will be cut by a machine word.
           A new initial state 0 is used (provided it is not already there).
           map is a memory area of m integers where the mapping of states is
           returned, if the caller wants to use it for something else.
        */

void regularRemapStates (Mask *trans, int m, Mask initial, Mask final,
                  Mask active, Mask *B, Mask **rtrans, Mask *rB, int *rm,
                  int *width, int *slices, Mask *rinitial, Mask *rfinal,
                  int **Map);

        /* reverses all the arrows of the NFA, exchanges initial and
           final states. the eps closures must have been done already */

void regularReverseArrows (Mask *trans, int m, Mask initial, 
                  Mask final, Mask **rtrans, Mask *rinitial, Mask *rfinal);

        /* builds a deterministic table for trans. it can handle multiwords
           in the det table too, but this is compatible with monoword usage.
           it receives the width of the deterministic table. the NFA is
           built so that the slices are always inside a single word */

void regularMakeDet (int width, Mask *trans, int m, Mask ***dtrans);

        /* same as regularMakeDet but the masks are simple */

void regularMakeDet1 (int width, Mask *trans, int m, mask ***dtrans);


        /* receives pat and its strlen m, the number of bits L to store the
           regexp, trees e and pos[]. It fills the preallocated tables B,S,A,
           and the automaton: transitions *trans (not preallocated), and
           initial and final states (preallocated) */

void regularLoadMasks (byte *pat, int m, int L, Tree *e, Tree **pos,
                       Mask *B, Mask *S, Mask A,
                       Mask **trans, Mask initial, Mask final);

        /* Finds a best factor to search the regexp, puts in *wlen the window 
	   length (0 if it recommends forward scanning). Returns the avg number
	   of chars inspected. In xm,xactive,xinitial,xfinal returns the
	   states selected for fast scanning, fwd and bwd verification. K is
	   a minimum number of chars that is read at each window */

double regularFindBest (Tree *e, Mask *trans, Mask *rtrans, int m,
                        Mask *B, Mask initial, Mask final, int *wlen,
                        Mask *dactive, Mask *dinitial, Mask *dfinal, int K);

        /* loads the fwd or bwd masks */                       

void regularLoadVerif (int m, Mask *B, int base, int sign, Mask **BV);

                /* loads masks for fast bwd/fwd */

regularScanData *regularLoadFast (int m, Mask *trans, int wlen, Mask *B,
                              Mask active, Mask initial, Mask final, int **map);
 
		/* frees */

void regularFreeScan (regularScanData *scan);

	/* the scanProc for regular expressions */
         
bool regularScan (byte **beg, byte **end, checkProc checkMatch, void *P,
		  regularScanData *scan);

#endif
