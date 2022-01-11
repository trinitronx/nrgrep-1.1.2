
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
 
#include "options.h"

        /* Global search options */
 
bool OptCaseInsensitive = false;   /* case insensitive search? */
byte *OptRecPatt = "\n";  /* record delimiter */
int OptRecChar = '\n';  /* OptRecPatt[0] if single char, else -1 */
bool OptRecPos = 1;  /* record starts at this position inside rec delimiter */
bool OptRecPrint = true;  /* print matching records? */
bool OptRecPrintFiles = false;  /* print whole matching files? */
bool OptRecFiles = false;  /* just report filenames? */
bool OptRecFileNames = true;  /* show filenames? */                             
int OptBufSize = 64*1024;  /* buffer size */
bool OptRecPositive = true;  /* show matching (not nonmatching) records */ 
bool OptRecNumber = false;  /* output record preceded by record number */
byte *OptRecSep = "";  /* output record separator */
bool OptLiteral = false;	/* take pattern literally? */
int OptErrors = 0;	/* number of errors permitted */
bool OptIns = true;         /* permit insertions in the pattern */
bool OptDel = true;         /* permit deletions in the pattern */
bool OptSubs = true;        /* permit substitutions in the pattern */
bool OptTransp = true;      /* permit transpositions in the pattern */          
bool OptWholeWord = false;   /* whole word matching */
bool OptWholeRecord = false;  /* occurrence has to match whole record */
bool OptStartLine = false;    /* occurrence has to start a line */
bool OptEndLine = false;    /* occurrence has to end a line */ 
int OptDetWidth = 16;	     /* width of deterministic tables
			        must be between 1 and W and divide W */
