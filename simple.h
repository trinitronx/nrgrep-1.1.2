
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

#ifndef SIMPLEINCLUDED
#define SIMPLEINCLUDED

	/* Search a simple pattern in a text */

#include "basics.h"	/* basics */
#include "options.h"	/* buffer management */
#include "parser.h"	/* parser */

typedef struct
   { mask B[256];	/* mask of forward or backward scan */
     int m;		/* length of subpattern used for scaning */
     bool bwd;		/* is it backward scan, right? */
   } simpleScanData;
   
typedef struct 
   { scanProc scanText; /* procedure for fast scanning */
     freeScanProc scanFree; /* procedure for freeing scanText */
     void *scanData;    /* data for fast scanning */
     int m;		/* length of subpattern used for scaning */
     Mask BfwdV[256];	/* mask for forward verification of the rest */
     Mask BbwdV[256];	/* mask for backward verification of the rest */
     int mbwdV;		/* length of extra backward verification */
     int mfwdV;		/* length of extra forward verification */
     bool raw;		/* used for raw or complete search */
   } simpleData;

	/* Preprocesses pat and creates a simpleData structure
	   for searching */

simpleData *simplePreproc (byte *pat, Tree *tree, Tree **pos);

	/* Frees P */

void simpleFree (simpleData *P);

        /* receives L in zero and returns there the number of bits to hold
           the pattern */

void simpleLength (Tree *tree, int *L);


	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

bool simpleSearch (byte **beg, byte **end, simpleData *P);

	/* Searches from *end-1 to *beg for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the last occurrence of P (first found in bwd scan) */
	/* This is used only to search the record separator */

bool simpleRevSearch (byte **beg, byte **end, simpleData *P);


        /* the scanProc for simple patterns. scanData is assumed to be
           of type simpleScan */

bool simpleScan (byte **beg, byte **end, checkProc checkMatch, void *P,
		 simpleScanData *scanData);

	/* loads mask B from pattern */

void simpleLoadMasks (byte *pat, int m, int L, Mask *B, Tree **pos);
	 
	/* Finds a best factor to search from *beg to *end-1, and
           recommends in *bwd a bwd o fwd search of it. Assumes that
           K extra characters are read in any window.  Returns avg
           work per text character. */                                         

double simpleFindBest (Mask *B, int m, int K,
		       bool *bwd, int *beg, int *end);
 
	/* loads bwd or fwd verification masks */

void simpleLoadVerif (int m, Mask *B, int base, int sign, Mask *BV);
	 
                /* loads masks for fast bwd/fwd */

simpleScanData *simpleLoadFast (bool bwd, Mask *B, int beg, int end);

                /* frees */

void simpleFreeScan (simpleScanData *scan);

#endif
