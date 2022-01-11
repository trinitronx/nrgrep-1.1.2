
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

#ifndef EEXTENDEDINCLUDED
#define EEXTENDEDINCLUDED

	/* Search an extended pattern in a text allowing errors */

#include "basics.h"	/* basics */
#include "options.h"	/* buffer management */
#include "extended.h"

typedef struct
   { extendedScanData *edata;	/* data for normal extended scanning */
     mask B0[256];	/* mask for exact search, initial states zeroed */
     int type;		/* type of search used: part in k+1, fwd or bwd scan */
     int m;		/* length of subpattern used for scaning */
     int k;		/* amount of errors */
     mask *ends;	/* ending positions for exact search */
     mask *V3;		/* some preallocated space */
   } eextendedScanData;

typedef struct 

   { escanProc scanText; /* procedure for fast scanning */
     freeScanProc scanFree; /* procedure for freeing scanText */
     void *scanData;    /* data for fast scanning */
     int k;		/* amount of errors */
     int type;		/* type of search used: part in k+1, fwd or bwd scan */
     Mask **BfwdV;	/* mask for forward verification of the rest */
     Mask **SfwdV;     /* S mask for forward verification of the rest */
     Mask *IfwdV,*FfwdV,*AfwdV; /* I,F,A masks for fwd verif of the rest */
     Mask *DifwdV;      /* initial D mask for forward verif of the rest */
     Mask **BbwdV;	/* mask for backward verification of the rest */
     Mask **SbwdV;     /* S mask for backward verification of the rest */
     Mask *IbwdV,*FbwdV,*AbwdV; /* I,F,A masks for bwd verif of the rest */
     Mask *DibwdV;      /* initial D mask for backward verif of the rest */
     int *mbwdV;	/* length of extra backward verification */
     int *mfwdV;	/* length of extra forward verification */
     Mask *V1;		/* some preallocated space */
     Mask V2;		/* some preallocated space */
   } eextendedData;

	/* Preprocesses pat and creates an eextendedData structure
	   for searching pat with k errors */

eextendedData *eextendedPreproc (byte *pat, Tree *tree, Tree **pos, int k);

	/* Frees P */

void eextendedFree (eextendedData *P);

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

bool eextendedSearch (byte **beg, byte **end, eextendedData *P);

	/* loads masks for fwd/bwd verification */

eextendedScanData *eextendedLoadFast (int wlen, int k, int type,
          Mask *B, Mask *S, Mask A, int *beg, int *end);

	/* frees */

void eextendedFreeScan (eextendedScanData *scan);

        /* the scanProc for eextended patterns. scanData is assumed to be
           of type eextendedScan */

bool eextendedScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
                    eextendedScanData *scan);

#endif
