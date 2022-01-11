
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

#ifndef EXTENDEDINCLUDED
#define EXTENDEDINCLUDED

	/* Search an extended pattern in a text */

#include "basics.h"	/* basics */
#include "options.h"	/* buffer management */
#include "parser.h"	/* parser */

typedef struct
   { mask B[256];	/* mask of forward or backward scan */
     mask S[256];	/* mask of staying states (for []*) */
     mask I,F,A;	/* mask of initial/final/? states (for []?) */
     int wlen;		/* window length, 0 => fwd scan */
     int m;		/* bitmask length */
   } extendedScanData;

typedef struct 
   { scanProc scanText; /* procedure for fast scanning */
     freeScanProc scanFree; /* procedure for freeing scanText */
     void *scanData;	/* data for fast scanning */
     mask *DV;		/* D mask preallocated for verification */
     Mask BfwdV[256];	/* mask for forward verification of the rest */
     Mask SfwdV[256];	/* S mask for forward verification of the rest */
     Mask IfwdV,FfwdV,AfwdV; /* I,F,A masks for forward verif of the rest */
     Mask DifwdV;	/* initial D mask for forward verif of the rest */
     Mask BbwdV[256];	/* mask for backward verification of the rest */
     Mask SbwdV[256];	/* S mask for backward verification of the rest */
     Mask IbwdV,FbwdV,AbwdV; /* I,F,A masks for backward verif of the rest */
     Mask DibwdV;	/* initial D mask for backward verif of the rest */
     int mbwdV;		/* length of extra backward verification */
     int mfwdV;		/* length of extra forward verification */
     int type;		/* FWD or BWD */
   } extendedData;

	/* Preprocesses pat and creates an extendedData structure
	   for searching */

extendedData *extendedPreproc (byte *pat, Tree *tree, Tree **pos);

	/* Frees P */

void extendedFree (extendedData *P);

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

bool extendedSearch (byte **beg, byte **end, extendedData *P);

        /* receives L in zero and returns there the number of bits to hold
           the pattern */

void extendedLength (Tree *tree, int *L);

        /* load masks B,S,A from pattern pat */

void extendedLoadMasks (byte *pat, int m, int L, Tree *tree, Tree **pos,
                        Mask *B, Mask *S, Mask *A);

        /* Finds a best factor to search from *beg to *end-1, and
           puts in *wlen the window length (0 if it recommends
           forward scanning). K is a number of chars always assumed
           to be scanned. Returns the avg number of chars inspected */         

double extendedFindBest (Mask *B, Mask *S, Mask A, int m, int K,
		         int *wlen, int *beg, int *end);
 
        /* loads the fwd or bwd masks */                       

void extendedLoadVerif (int m, Mask *B, Mask *S, Mask A, int base, int sign,
		       Mask *BV, Mask *SV, Mask *IV, Mask *FV, Mask *AV,
		       Mask *DiV);

                /* loads masks for fast bwd/fwd */

extendedScanData *extendedLoadFast (int wlen, Mask *B, Mask *S, Mask A, 
				    int beg, int end);

                /* frees */

void extendedFreeScan (extendedScanData *scan);

	/* the scanProc for extended patterns. scanData is assumed to be
	   of type extendedScan */

bool extendedScan (byte **beg, byte **end, checkProc checkMatch, void *P,
                   extendedScanData *scan);



#endif
