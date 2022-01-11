
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

#ifndef BUFFERINCLUDED
#define BUFFERINCLUDED

#include "basics.h"

	/* Buffer management */

typedef struct sBuffer
  { int size;		/* size of this buffer */
    int file;		/* associated file */
    int fpos;		/* pos of current buffer in file */
    byte *data;		/* current text */
    int dsize;		/* size of current text */
  } *Buffer;

	/* Creates a new buffer */

Buffer bufCreate (void);

	/* Terminates use of buffer B */

void bufDestroy (Buffer B);

	/* Assigns a file to B */

void bufSetFile (Buffer B, int file);

        /* Gets buffer size of B */                                             

int bufSizer (Buffer B);
 
	/* Reads a new buffer starting at position pnext */

void bufLoad (Buffer B, byte *pnext);

  	/* Tells whether there is no more data in B + its file */

bool bufEmpty (Buffer B);

  	/* Tells whether there is no more data to read in B's file */

bool bufEof (Buffer B);

	/* Tells the file position of the current buffer */

int bufTextPos (Buffer B);

	/* Gives the memory area of the current buffer */
	/* One text position before the given buffer can be touched */

void bufCurrent (Buffer B, byte **base, byte **top);

#endif
