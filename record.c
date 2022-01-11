
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

	/* Drives the search considering the records */

	/* A match totally overlapping a separator is disregarded. */

	/* A record longer than the buffer size can be cut at arbitrary
	   points and some matches can be missed. */

#include "record.h"
#include "options.h"

	/* module variables */

static simpleData *RecPatt = NULL;
static int RecLen = 0;
static bool RecByRec,RecBeg,RecEnd;

void recPreproc (void)

		/* makes the preprocessing for record handling */
		 
   { Tree *tree;
     Tree **pos;
     int m = strlen(OptRecPatt);
     int i;
     int tmpOptRecChar;
     pos = malloc (m * sizeof(Tree*));
     tree = parse (OptRecPatt, m, pos);
     if ((tree == NULL) || tree->eps || (detClass(tree,NULL,0) != SIMPLE))
	error1("Syntax error in delimiter %s",OptRecPatt);
     tmpOptRecChar = OptRecChar;
     OptRecChar = -1;
     RecPatt = simplePreproc (OptRecPatt,tree,pos);
     OptRecChar = tmpOptRecChar;
     RecPatt->raw = true;
     RecLen = 0;
     simpleLength(tree,&RecLen);
     freeTree (tree); free (pos);
   }

void recFree (void)

		/* frees the structures of record handling */
		 
   { if (RecPatt) simpleFree (RecPatt);
     RecPatt = NULL; RecLen = 0;
   }

void recGetRecord (byte *ptr, byte *bbeg, byte *bend,
		   byte **rbeg, byte **rend, byte **obeg, byte **oend)

		/* ptr is a point inside a buffer [bbeg,bbend-1]. 
		   this returns [*rbeg,*rend-1] as the record containing ptr. 
		   it also returns [*obeg,*oend-1] as the record to return
		   if the match is successful. */

   { byte *beg,*end;
     if (RecByRec)   /* buffer is already the record */
	{ *rbeg = bbeg;
	  *rend = bend;
  	  *obeg = bbeg + (RecBeg ? - RecLen + OptRecPos : 0);
	  *oend = bend + (RecEnd ? OptRecPos : 0);
	  return;
	}
     
     	   /* first, reverse search of record beginning */

     beg = bbeg; end = ptr + RecLen - 1;
     if (end > bend) end = bend;
     if (simpleRevSearch (&beg,&end,RecPatt))
        { *rbeg = end;
	  *obeg = end - RecLen + OptRecPos;
        }
     else *rbeg = *obeg = bbeg;

     	   /* second, forward search of record ending */

     beg = ptr; end = bend;
     if (simpleSearch (&beg,&end,RecPatt))
        { *rend = beg;
	  *oend = beg + OptRecPos;
        }
     else *rend = *oend = bend;
   }

	        /* check that a pattern occurrence in [pbeg,pend-1]
	           is legal inside the record [rbeg,rbeg-1] because of
		   context restrictions */

bool recCheckLeftContext (byte *pbeg, byte *rbeg)

   { if (OptWholeWord)  /* match whole words */
	{ if ((pbeg > rbeg) && (isalnum(pbeg[-1]))) return false;
	}
     if (OptStartLine)  /* occurrence must start a line */
	{ if ((pbeg > rbeg) && (pbeg[-1] != '\n')) return false;
	}
     if (OptWholeRecord)  /* occurrence has to be the whole record */
	{ if (pbeg > rbeg) return false;
	}
     return true;
   }

bool recCheckRightContext (byte *pend, byte *rend)

   { if (OptWholeWord)  /* match whole words */
	{ if ((pend < rend) && (isalnum(pend[0]))) return false;
	}
     if (OptEndLine)  /* occurrence must end a line */
	{ if ((pend < rend) && (pend[0] != '\n')) return false;
	}
     if (OptWholeRecord)  /* occurrence has to be the whole record */
	{ if (pend < rend) return false;
	}
     return true;
   }

int recSearchFile (char *fname, Buffer B, searchData *P)

        /* searches the file handled by B for P using R as record
	   separator, reports matches as appropriate and returns number
	   of matches */
 
   { byte *text,*top,*nextbuf,*pbeg,*pend;
     int matches = 0;
     bool firstReport = true;

     RecByRec = false;

     while (!bufEmpty(B))
	{ bufCurrent (B,&text,&top);
		/* find last record limit */
	  if (bufEof(B)) nextbuf = top; /* end of record = end of file */
	  else 
	     { pbeg = text; pend = top;
	       if (simpleRevSearch (&pbeg,&pend,RecPatt) && (pbeg != text)) 
	          { nextbuf = pbeg;
	            top = pend;
	          }
	       else
	          {    /* record will be split, warn */
	            warn1 ("Record longer than buffer size (%i) has been split",
		           bufSize(B));
	            nextbuf = top;
	          }
	     }
		/* now search the pattern in the buffer */
	  while (true)
	     { pbeg = text; pend = top;
	       if (searchScan (&pbeg,&pend,P))
		  {   /* report the matching record */
		    if (OptRecFiles)  /* just report the file */
		       { if (OptRecFileNames && (fname != NULL))
			    printf ("%s\n",fname);
			 return 1;
		       }
		    if (OptRecFileNames && (fname != NULL) && firstReport)
		       { printf ("%s:",fname);
			 if (OptRecPrint) printf ("\n");
			 firstReport = false;
		       }
		    if (OptRecPrint) /* show the matches */
		       { fwrite (pbeg,pend-pbeg,1,stdout);
			 printf (OptRecSep);
		       }
		    matches++; /* count the matches */
		    if (pend == top) break;  /* buffer totally processed */
		    text = pend - OptRecPos + RecLen;
		  }
	       else break;  /* no more occurrences in this buffer */
	     }
		/* load next buffer */
	  bufLoad (B,nextbuf);
	}
     if (OptRecFileNames && (fname != NULL) && !OptRecPrint && (matches > 0)) 
	printf(" %i\n",matches);
     return matches;
   }

int recSearchRecFile (char *fname, Buffer B, searchData *P)

        /* searches the file handled by B for P using R as record
	   separator, but scans record by record. this is our way to
	   handle OptRecNumber or !OptRecPositive */
 
   { byte *text,*top,*rbeg,*rend,*pbeg,*pend,*recb,*rece;
     byte *lastrep,*rbeg2,*rend2;
     int matches = 0;
     bool firstReport = true;
     int recNum = 0;

     RecByRec = true;

     while (!bufEmpty(B))
	{ bufCurrent (B,&text,&top);
		/* go record by record */
	  rbeg = text-1;
	  while (true)
	     { recb = ++rbeg; rend = top;
	       if (!simpleSearch (&rbeg,&rend,RecPatt))
	          { if (recb <= text+1)
	               { if (!bufEof(B)) /* record will be split, warn */
			    warn1 ("Record longer than buffer size (%i)"
				   " has been split", bufSize(B));
			 rbeg = top;
		       }
		    else break;   /* read another buffer */
		  }
	       rece = rbeg;
	       if (rbeg == text) continue;/* rec delim at the beginning, skip */
	       if (recb != text) recb += RecLen-1;
			/* now we search the pattern inside the record */
	       recNum++;
	       RecBeg = (recb != text);
	       RecEnd = (rece != top);
	       pbeg = recb; pend = rece;
	       if (searchScan (&pbeg,&pend,P) == OptRecPositive)
		  {   /* found it, report as appropriate */
		    if (OptRecFileNames && (fname != NULL) && firstReport)
		       { printf ("%s:",fname);
			 if (OptRecPrint) printf ("\n");
			 firstReport = false;
		       }
		    if (OptRecPrintFiles) return 1;
		    if (OptRecNumber) printf ("%i:%s",recNum,OptRecSep);
		    if (OptRecPrint)
	               { pbeg = recb + (RecBeg ? - RecLen + OptRecPos : 0);
	                 pend = rece + (RecEnd ? OptRecPos : 0);
		         fwrite (pbeg,pend-pbeg,1,stdout);
		         printf (OptRecSep);
		       }
		    matches++;
		  }
	     }
	    /* load next buffer */
	   bufLoad (B,recb-1);
	}
     if (OptRecFileNames && (fname != NULL) && !OptRecPrint && (matches > 0)) 
	printf(" %i\n",matches);
     return matches;
     return matches;
   }

