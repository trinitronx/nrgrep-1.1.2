
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

#ifndef EREGULARINCLUDED
#define EREGULARINCLUDED

	/* Search a regular expression in a text allowing errors */

#include "basics.h"	/* basics */
#include "options.h"	/* buffer management */
#include "regular.h"

typedef struct
   { regularScanData *edata;  /* data for normal regular scanning */
     int type;		/* type of search used: part in k+1, fwd or bwd scan */
     int k;		/* number of errors */
     int dm;		/* with of final masks */
     mask *V3;		 /* extra space for searching */
   } eregularScanData;

typedef struct 

   { escanProc scanText; /* procedure for fast scanning */
     freeScanProc scanFree; /* procedure for freeing scanText */
     void *scanData;    /* data for fast scanning */
     int type;		/* type of search used: part in k+1, fwd or bwd scan */
     int k;		/* number of errors */
     int m,dm;          /* number of states orig and scanning autom */

     mask match;	/* matched in scanning autom */
     int width;		/* width of slices for verif autom */
     Mask **dftransV;   /* deterministic version of fwd verif automaton */
     Mask **dbtransV;   /* deterministic version of bwd verif automaton */
     Mask initial,final; /* initial,final states for verif */
     Mask B[256]; 	 /* B masks for verification */
     int *revMap;	 /* map local -> global states */
     Mask V1;		 /* extra space for verification */
     Mask *V2;		 /* extra space for verification */
   } eregularData;

	/* Preprocesses pat and creates an eregularData structure
	   for searching */

eregularData *eregularPreproc (byte *pat, Tree *tree, Tree **pos, int k);

	/* Frees P */

void eregularFree (eregularData *P);

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

bool eregularSearch (byte **beg, byte **end, eregularData *P);

                /* loads masks for fast bwd/fwd */

eregularScanData *eregularLoadFast (int m, Mask *trans, int wlen,
                                   Mask *B, Mask *active, Mask *initial,
                                   Mask *final, int k, int type, int **map);

		/* frees */

void eregularFreeScan (eregularScanData *scan);

	/* the scanProc for eregular expressions */
         
bool eregularScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		  eregularScanData *scan);

#endif
