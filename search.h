
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

#ifndef SEARCHINCLUDED
#define SEARCHINCLUDED

	/* Search a pattern in a text, general interface */

#include "basics.h"	/* basics */
#include "options.h"	/* buffer management */

typedef struct 

   { int searchType;	/* type of pattern */
     void *preprocData;	/* data from preprocessing */
     bool (*search)();  /* search function */
   } searchData;

	/* Preprocesses pat and creates a searchData structure
	   for searching */

searchData *searchPreproc (byte *pat);

	/* Frees P */

void searchFree (searchData *P);

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the printable record containing the occurrence of P */

bool searchScan (byte **beg, byte **end, searchData *P);

#endif
