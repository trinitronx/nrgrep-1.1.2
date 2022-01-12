
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

#ifndef BASICSINCLUDED
#define BASICSINCLUDED

	/* some basic types */

typedef unsigned char byte;

typedef int bool;
#define true 1
#define false 0

	/* some basic operations */

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

	/* some public modules */

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>  /* for exit, malloc, realloc, free */
#include <string.h>  /* for strlen */

        /* some constants */

#define KPLUS1 1
#define BWD    2
#define FWD    3

	/* some parameterized procedures */

        /* Checks from *beg to *end-1 whether there is a match around pos and 
	   returns in *beg and *end the limits. "Around" means that pos is the 
	   initial position of the search pattern given to scanProc.
        */

typedef bool (*checkProc)(void *P,byte *pos,byte **beg, byte **end);
typedef bool (*echeckProc)(void *P,int i,byte *pos,byte **beg, byte **end);

        /* Searches from *beg to *end-1 for P. Returns if it could find
           it. In case of returning true, *beg and *end are set to limit
           the first occurrence of P */

typedef bool (*scanProc)(byte **beg, byte **end,checkProc checkMatch,
			 void *P,void *scanData);
typedef bool (*escanProc)(byte **beg, byte **end,echeckProc checkMatch,
			  void *P,void *scanData);

	/* Frees scanData */

typedef void (*freeScanProc)(void *scanData);

	/* some private modules */

#include "except.h"	/* exception handling */
#include "memio.h"	/* memory and i/o management */
#include "bitmasks.h"	/* bitmask management */

#endif
