
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

#ifndef RECORDINCLUDED
#define RECORDINCLUDED

	/* Search the record separator */

	/* A match totally containing a separator is disregarded. */

	/* A delimiter overlapping a match is disregarded */

#include "basics.h"	/* basics */
#include "buffer.h"	/* buffer management */
#include "simple.h"	/* simple searching */
#include "search.h"	/* general searching */
 
                /* makes the preprocessing for record handling */

void recPreproc (void);
	 
	        /* frees the structures of record handling */

void recFree (void);
	 
		/* ptr is a point inside a buffer [bbeg,bbend-1].
		   this returns [*rbeg,*rend-1] as the record containing ptr.
		   it also returns [*obeg,*oend-1] as the record to return
		   if the match is successful. */

void recGetRecord (byte *ptr, byte *bbeg, byte *bbend,
		   byte **rbeg, byte **rend, byte **obeg, byte **oend);

	        /* check that a pattern occurrence in [pbeg,pend-1]
	           is legal inside the record [rbeg,rbeg-1] because of
	           context restrictions */

bool recCheckLeftContext (byte *pbeg, byte *rbeg);
bool recCheckRightContext (byte *pend, byte *rend);
	 

        /* searches the file handled by B for P using R as record
           separator, reports matches as appropriate and returns number
           of matches */                                                        

int recSearchFile (char *fname, Buffer B, searchData *P);

        /* searches the file handled by B for P using R as record
           separator, but scans record by record. this is our way to
           handle OptRecNumber or !OptRecPositive */                               

int recSearchRecFile (char *fname, Buffer B, searchData *P);

#endif
