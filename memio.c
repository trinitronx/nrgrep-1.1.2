
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

/* #include "memio.h" cannot be done because of circular redefinition */

#include <stdio.h>
#include "except.h"

void *mymalloc (int n)

	/* a more robust malloc routine */

  { void *p;
    if (n == 0) return NULL;
    p = (void *) malloc (n);
    if (p == NULL) error1("cannot allocate %i bytes",n);
    return p;
  }

void myfree (void *p)

	/* a more robust free routine */

  { if (p != NULL) free (p);
  }

void *myrealloc (void *p, int n)

	/* a more robust realloc routine */

  { if (n == 0)
       { myfree (p);
	 return NULL;
       }
    if (p == NULL) return mymalloc (n);
    p = (void*) realloc (p,n);
    if (p == NULL) error1("cannot reallocate to %i bytes",n);
    return p;
  }

int myread (int file, void *buf, int count)

	/* a more robust read routine */

  { int n;
    if (count == 0) return 0;
    n = read (file,buf,count);
    if (n == -1) error1("cannot read %i bytes from file",count);
    return n;
  }

int mywrite (int file, void *buf, int count)

	/* a more robust write routine */

  { int n;
    if (count == 0) return 0;
    n = write (file,buf,count);
    if (n == -1) error1("cannot write %i bytes to file",count);
    return n;
  }

