
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

#ifndef EXCEPTINCLUDED
#define EXCEPTINCLUDED

	/* Exception handling module */

#include <errno.h>

#define error0(str) \
    { fprintf(stderr,"Fatal error: "); \
      fprintf(stderr,str); \
      fprintf (stderr,"\n -- errno %i, line %i of %s\n\n", \
	       errno,__LINE__,__FILE__); \
      exit(1); \
    }

#define error1(str,p1) \
    { fprintf(stderr,"Fatal error: "); \
      fprintf(stderr,str,p1); \
      fprintf (stderr,"\n -- errno %i, line %i of %s\n\n", \
	       errno,__LINE__,__FILE__); \
      exit(1); \
    }

#define error2(str,p1,p2) \
    { fprintf(stderr,"Fatal error: "); \
      fprintf(stderr,str,p1,p2); \
      fprintf (stderr,"\n -- errno %i, line %i of %s\n\n", \
	       errno,__LINE__,__FILE__); \
      exit(1); \
    }

#define error3(str,p1,p2,p3) \
    { fprintf(stderr,"Fatal error: "); \
      fprintf(stderr,str,p1,p2,p3); \
      fprintf (stderr,"\n -- errno %i, line %i of %s\n\n", \
	       errno,__LINE__,__FILE__); \
      exit(1); \
    }

#define warn0(str) \
    { fprintf(stderr,"Warning: "); \
      fprintf(stderr,str); \
      fprintf(stderr,"\n"); \
    }

#define warn1(str,p1) \
    { fprintf(stderr,"Warning: "); \
      fprintf(stderr,str,p1); \
      fprintf(stderr,"\n"); \
    }

#define warn2(str,p1,p2) \
    { fprintf(stderr,"Warning: "); \
      fprintf(stderr,str,p1,p2); \
      fprintf(stderr,"\n"); \
    }

#define warn3(str,p1,p2,p3) \
    { fprintf(stderr,"Warning: "); \
      fprintf(stderr,str,p1,p2,p3); \
      fprintf(stderr,"\n"); \
    }

#endif
