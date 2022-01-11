
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
 
#ifndef PARSERINCLUDED
#define PARSERINCLUDED

#include "basics.h"
 
        /* Gets one character assuming it is not special.
           Returns -1 in case of syntax error: \x not followed by 2 hex
           digits or \ terminating the string. */
        /* Syntax: \n -> newline, \t -> tab, \xXY -> hex ascii code XY,
                   \X -> literal character X, else normal char */

int getAchar (byte *pat, int *i);

        /* Gets a class of characters from pat[i...] (sets i at the next
           position to read) and sets the appropriate bits in B[*] column l. 
	   From outside it should be used to work on positions marked by
	   the pos[] array returned by parse */

bool getAclass (byte *pat, int *i, Mask *B, int l);

	/* probabilities of individual letters */

extern double letterProb[256];

#define STR 0
#define STAR 1
#define OOR 2
#define CONC 3
#define QUESTION 4
#define PLUS 5

typedef struct sTree
   { int type;		   /* node type: STR,STAR,OOR,CONC,QUESTION,PLUS */
     int pos;		   /* for STR type, position in bit mask */
     bool eps;		   /* subexpression recognizes epsilon? */
     struct sTree *e1,*e2; /* subexpressions */
     Mask maskPos;	   /* positions covered by subexpression */
     Mask firstPos;	   /* initial positions of subexpression */
     Mask lastPos;	   /* final positions of subexpression */
   } Tree;

        /* parse makes a pre parsing of the regular expression. It
	   receives the pattern, its strlen m, and an array pos of size m.
	   It parses pat and returns NULL in case of
	   parsing error. Otherwise it returns a syntax tree whose root tells
	   the type of pattern, and in pos[i] it puts a pointer to the node
	   where the class starting at pat[i] must be set.
	*/
           
Tree *parse (byte *pat, int m, Tree **pos);

        /* determines the class of pattern for a subset active (all if NULL) */

	/* detClass only returns types 1,2,3. order is important */

#define SIMPLE 1
#define EXTENDED 2
#define REGULAR 3
#define ESIMPLE 4
#define EEXTENDED 5
#define EREGULAR 6

int detClass (Tree *e, Mask active, int m);

        /* compute maskPos (positions of subexpression) */

void setMaskPos (Tree *e, int L);

        /* frees a structure allocated by parse */

void freeTree (Tree *e);

#endif
