
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

#ifndef ESIMPLEINCLUDED
#define ESIMPLEINCLUDED

	/* Search a simple pattern in a text */

#include "basics.h"	/* basics */
#include "options.h"	/* buffer management */
#include "simple.h"	/* search w/o errors */

typedef struct
   { simpleScanData *sdata;  /* the data for a simple scan */
     mask B0[256];	/* mask for exact search, initial states zeroed */
     int k;		/* amount of errors */
     int type;		/* type of search used: part in k+1, fwd or bwd scan */
     mask *V3;		/* some preallocated space */
   } esimpleScanData;

typedef struct 

   { escanProc scanText; /* procedure for fast scanning */
     freeScanProc scanFree; /* procedure for freeing scanText */
     void *scanData;    /* data for fast scanning */
     int k;		/* amount of errors */
     int type;		/* type of search used: part in k+1, fwd or bwd scan */
     Mask **BfwdV;	/* mask for forward verification of the rest */
     Mask **BbwdV;	/* mask for backward verification of the rest */
     int *mbwdV;	/* length of extra backward verification */
     int *mfwdV;	/* length of extra forward verification */
     Mask *V1;		/* some preallocated space */
     Mask V2;		/* some preallocated space */
   } esimpleData;

	/* Preprocesses pat and creates an esimpleData structure
	   for searching pat with k errors */

esimpleData *esimplePreproc (byte *pat, Tree *tree, Tree **pos, int k);

	/* Frees P */

void esimpleFree (esimpleData *P);

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

bool esimpleSearch (byte **beg, byte **end, esimpleData *P);

	/* loads masks for fast fwd/bwd scanning */

esimpleScanData *esimpleLoadFast (int wlen, int k, int type, Mask *B, int *beg);

	/* frees */

void esimpleFreeScan (esimpleScanData *scan);

        /* the scanProc for esimple patterns. scanData is assumed to be
           of type esimpleScan */

bool esimpleScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
                  esimpleScanData *scan);

#endif
