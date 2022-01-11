
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

#ifndef OPTIONSINCLUDED
#define OPTIONSINCLUDED

#include "basics.h"

	/* Global search options */

extern bool OptCaseInsensitive;   /* case insensitive search? */
extern byte *OptRecPatt;  /* record delimiter */
extern int OptRecChar;
extern bool OptRecPos; /* record starts at this position inside rec delimiter */
extern bool OptRecPrint;  /* print matching records? */
extern bool OptRecPrintFiles;  /* print whole matching files? */
extern bool OptRecFiles;  /* just report filenames? */
extern bool OptRecFileNames;  /* show filenames? */
extern int OptBufSize;  /* buffer size */ 
extern bool OptRecPositive;  /* show matching (not nonmatching) records */
extern bool OptRecNumber;  /* output record preceded by record number */
extern byte *OptRecSep;  /* output record separator */
extern bool OptLiteral;	  /* take the pattern literally */
extern int OptErrors;      /* number of errors permitted */
extern bool OptIns;	    /* permit insertions in the pattern */
extern bool OptDel;	    /* permit deletions in the pattern */
extern bool OptSubs;	    /* permit substitutions in the pattern */
extern bool OptTransp;	    /* permit transpositions in the pattern */
extern bool OptWholeWord;   /* whole word matching */
extern bool OptWholeRecord;  /* occurrence has to match whole record */
extern bool OptStartLine;    /* occurrence has to start a line */
extern bool OptEndLine;    /* occurrence has to end a line */
extern int OptDetWidth;	    /* deterministic width for tables */

#endif
