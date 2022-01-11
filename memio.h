
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

#ifndef MEMIOINCLUDED
#define MEMIOINCLUDED

	/* Memory and I/O management module */

void *mymalloc (int n);
void myfree (void *p);
void *myrealloc (void *p, int n);

#define malloc(n) mymalloc(n)
#define free(p) myfree(p)
#define realloc(p,n) myrealloc(p,n)

int myread (int file, void *buf, int count);
int mywrite (int file, void *buf, int count);

#define read(f,b,c) myread(f,b,c)
#define write(f,b,c) mywrite(f,b,c)

#endif
