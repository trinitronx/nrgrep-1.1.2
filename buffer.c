
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

	/* Buffer management */

#include "buffer.h"
#include "options.h"

Buffer bufCreate (void)

	/* Creates a new buffer */

  { Buffer B = malloc (sizeof(struct sBuffer));
    B->size = OptBufSize;
    B->data = malloc (1+B->size); B->data++;
    B->dsize = 0;
    return B;
  }

void bufDestroy (Buffer B)

	/* Terminates use of buffer B */

  { free (B->data-1);
    free (B);
  }

void bufSetFile (Buffer B, int file)

	/* Sets the file for buffer B */

  { B->file = file;
    B->fpos = 0;
    B->dsize = read (B->file,B->data,B->size);
  }

int bufSize (Buffer B)

	/* Gets buffer size of B */

  { return B->size;
  }

void bufLoad (Buffer B, byte *pnext)

	/* Reads a new buffer starting at position next */

  { register byte *src = pnext;
    register byte *dst = B->data;
    register int mov = B->dsize - (src - dst);

	/* move unused part to the beginning (memcpy is not portable) */
    B->dsize = mov;
    while (mov--) *dst++ = *src++;
    mov = B->size - B->dsize;
    B->dsize += read (B->file, dst, mov);
    B->fpos += mov;
  }

bool bufEmpty (Buffer B)

  	/* Tells whether there is no more data in B + its file */

  { return B->dsize == 0;
  }

bool bufEof (Buffer B)

  	/* Tells whether there is no more data to read in B's file */

  { return B->dsize < B->size;
  }

int bufTextPos (Buffer B)

	/* Tells the file position of the current buffer */

  { return B->fpos;
  }

void bufCurrent (Buffer B, byte **base, byte **top)

	/* Gives the memory area of the current buffer */

  { *base = B->data;
    *top = B->data + B->dsize;
  }

