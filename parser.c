
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

#include "parser.h"
#include "options.h"

#define Sep '#'
#define All '.'
#define OpenPar '('
#define ClosePar ')'
#define Star '*'
#define Plus '+'
#define Question '?'
#define Or '|'
#define OpenBracket '['
#define CloseBracket ']'
#define BackSlash '\\'
#define Dash '-'
#define Neg '^'

int getAchar (byte *pat, int *i)
 
        /* Gets one character assuming it is not special.
	   Returns -1 in case of syntax error: \x not followed by 2 hex
	   digits or \ terminating the string. */
	/* Syntax: \n -> newline, \t -> tab, \xXY -> hex ascii code XY,
		   \X -> literal character X, else normal char */
 
   { int c;
     if (pat[*i] == BackSlash)
        { (*i)++;
          if (pat[*i] == 'n') c = '\n';
          else if (pat[*i] == 't') c = '\t';
          else if (tolower(pat[*i]) == 'x')
                  { (*i)++;
                     if (isdigit(pat[*i])) c = pat[*i] - '0';
                     else if (isxdigit(pat[*i])) c = (tolower(pat[*i])-'a')+10;
		     else return -1;
                     (*i)++; c <<= 4;
                     if (isdigit(pat[*i])) c += pat[*i] - '0';
                     else if (isxdigit(pat[*i])) c += (tolower(pat[*i])-'a')+10;
		     else return -1;
                  }
	  else if (pat[*i]) c = pat[*i];
	  else return -1;
        }
     else c = pat[*i];
     (*i)++;
     return c;
   }

bool getAclass (byte *pat, int *i, Mask *B, int l)

	/* Gets a class of characters from pat[i...] (sets i at the next
	   position to read) and sets the appropriate bits in
	   B[*] column l. It assumes that the first character is not a
	   metacharacter. It returns true if it could read a class, false
	   in case of syntax error: from getAchar or missing ]. It can
	   also be used to test syntax, with B = NULL */
	/* Syntax: class -> [set*] or [^set*] (^ complements the set)
		   set -> X-Y  (range of chars) or X (single char)
		   X,Y -> . (a class with all), # (a separator) or getAchar
	*/

   { int c,c1,c2;
     bool pos = true;
     switch (pat[*i])
       { case OpenBracket: /* a class is opened */
            (*i)++;
            if (pat[*i] == Neg)  /* reverse meaning */
               { if (B) for (c=0; c<256; c++) SET(B[c],l);
                 pos = false;
                 (*i)++;
               }
            while (pat[*i] && (pat[*i] != CloseBracket))
               { c1 = getAchar (pat,i); if (c1 == -1) return false;
                 if ((pat[*i] == Dash) && pat[*i+1] && 
		     (pat[*i+1] != CloseBracket))
                    { (*i)++;
                      c2 = getAchar (pat,i); if (c2 == -1) return false;
                    }
                 else c2 = c1;
                 if (B)
		    for (c=c1; c<=c2; c++)
                       { if (pos) SET(B[c],l);
                         else CLEAR(B[c],l);
                       }
               }
	    if (!pat[*i]) return false;  /* ] missing */
            (*i)++;
            break;
         case All:       /* a class with everything */
            if (B) for (c=0; c<256; c++) SET(B[c],l);
            (*i)++;
            break;
         case Sep:       /* a class of separators */
            if (B) for (c=0; c<256; c++) if (!isalnum(c)) SET(B[c],l);
            (*i)++;
            break;
         case OpenPar: case ClosePar: 
         case Or: case Question: case Star: case Plus:
                        /* higher level constructions, reject */
            return false;
            break;
         default:  /* a plain letter or escape code */
            c = getAchar (pat,i); if (c == -1) return false;
            if (B) SET(B[c],l);
            break;
       }

            /* expand for case insensitive searching */
     if (B && OptCaseInsensitive)
        { for (c='A'; c<='Z'; c++)
              if (pos && ISSET(B[c],l)) SET(B[tolower(c)],l);
              else if (!pos && !ISSET(B[c],l)) CLEAR(B[tolower(c)],l);
          for (c='a'; c<='z'; c++)
              if (pos && ISSET(B[c],l)) SET(B[toupper(c)],l);
              else if (!pos && !ISSET(B[c],l)) CLEAR(B[toupper(c)],l);
        }

     return true;
   }

double letterProb[256] = {
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000344, 0.020793, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.146588, 0.000043, 0.000460, 0.000398, 0.011430, 0.003034, 0.001013, 0.001707, 
0.004156, 0.004162, 0.000506, 0.000998, 0.008441, 0.003342, 0.009616, 0.000903, 
0.002255, 0.004002, 0.002441, 0.001222, 0.000937, 0.001102, 0.000874, 0.000828, 
0.000970, 0.001810, 0.000679, 0.000168, 0.000190, 0.001562, 0.000143, 0.000035, 
0.000086, 0.002093, 0.001334, 0.001530, 0.000818, 0.000981, 0.001181, 0.000571, 
0.000754, 0.001534, 0.000156, 0.000228, 0.000656, 0.001308, 0.000922, 0.001299, 
0.001202, 0.000261, 0.000689, 0.001809, 0.003403, 0.000669, 0.000340, 0.000961, 
0.000158, 0.000390, 0.000234, 0.000847, 0.015840, 0.000846, 0.001258, 0.001695, 
0.000715, 0.053857, 0.011376, 0.027900, 0.021596, 0.094887, 0.015707, 0.013246, 
0.030408, 0.054368, 0.000933, 0.003729, 0.028211, 0.020693, 0.048064, 0.047054, 
0.018812, 0.002436, 0.044806, 0.048118, 0.065831, 0.016154, 0.006572, 0.008692, 
0.005656, 0.007099, 0.001124, 0.008146, 0.000445, 0.008146, 0.001852, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000035, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000019, 0.000000, 0.000000, 0.000000, 0.000032, 0.000000, 0.000000, 
0.000000, 0.000009, 0.000000, 0.000040, 0.000000, 0.000000, 0.000000, 0.000000, 
0.000000, 0.000000, 0.000027, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
};


	/* Simple syntax:
		inside a class [...]
			^ at the beginning: negation of class
			\c: take character c as literal
			\n,\t: newline, tab
			\xdd: hex ascii code
		outside a class, in addition
 			# = separator
			. = a class of all
	*/

	/* Extended Syntax:
		inside a class [...]
			^ at the beginning: negation of class
			\c: take character c as literal
			\n,\t: newline, tab
			\xdd: hex ascii code
	        outside a class, in addition
	                <class>+
	                <class>?
			<class>*
			# = separator
	                . = a class of all                                      
	*/

	/* Regexp Syntax:
		inside a class [...]
			^ at the beginning: negation of class
			\c: take character c as literal
			\n,\t: newline, tab
			\xdd: hex ascii code
	        outside a class, in addition
	                <class>+
	                <class>?
			<class>*
			# = separator
	                . = a class of all                                      
		unique of regexps
			(exp), exp|exp, exp+, exp*, exp?, exp exp
	*/


void freeTree (Tree *e)

	/* frees a structure allocated by parse */
 
  { if (e == NULL) return;
    free (e->maskPos); free (e->firstPos); free (e->lastPos);
    freeTree (e->e1); freeTree (e->e2);
    free (e);
  }

static Tree *parseOr (byte *pat, int *i, Tree **pos);

static Tree *parseOne (byte *pat, int *i, Tree **pos)

	/* parses one subexpression */

  { Tree *e;
    if (pat[*i] == OpenPar)  /* opens subexpression */
       { (*i)++;
	 e = parseOr (pat,i,pos);
	 if (e == NULL) return NULL;
	 if (pat[*i] != ClosePar) { freeTree (e); return NULL; }
	 (*i)++;
       }
    else 
       { e = malloc (sizeof(Tree));
	 e->type = STR;
	 e->e1 = e->e2 = NULL;
		/* could be the empty string if higher level chars are here */
	 e->eps = !pat[*i] || (pat[*i] == ClosePar) || (pat[*i] == Star) ||
		  (pat[*i] == Plus) || (pat[*i] == Question) || (pat[*i] == Or);
	 if (!e->eps)
	    { pos[*i] = e;
	      if (!getAclass(pat,i,NULL,0)) { free (e); return NULL; }
	    }
	 e->maskPos = e->firstPos = e->lastPos = NULL;
       }
    return e;
  }

static Tree *parseOneClosed (byte *pat, int *i, Tree **pos)

	/* parses one subexpression plus closures +,*,? */

  { Tree *e,*e1;

    e = parseOne (pat,i,pos); if (e == NULL) return NULL;

    while (true)
        switch (pat[*i])  
           { case Plus:   /* plus */
	        e1 = malloc (sizeof(Tree));
		e1->type = PLUS; e1->eps = e->eps;
		e1->e1 = e; e1->e2 = NULL;
		e = e1;
	        e->maskPos = e->firstPos = e->lastPos = NULL;
	        (*i)++;
	        break;
             case Star:   /* star */
	        e1 = malloc (sizeof(Tree));
	        e1->type = STAR; e1->eps = true;
	        e1->e1 = e; e1->e2 = NULL;
	        e = e1;
	        e->maskPos = e->firstPos = e->lastPos = NULL;
	        (*i)++;
	        break;
             case Question:   /* previous is optional */
		e1 = malloc (sizeof(Tree));
		e1->type = QUESTION; e1->eps = true;
		e1->e1 = e; e1->e2 = NULL;
		e = e1;
	        e->maskPos = e->firstPos = e->lastPos = NULL;
	        (*i)++;
	        break;
	     default: return e;  /* end of closed expression */
           }
      }

static Tree *parseConc (byte *pat, int *i, Tree **pos)

	/* parses a concatenation */

  { Tree *e,*e1;
    e = parseOneClosed (pat,i,pos); if (e == NULL) return NULL;
    switch (pat[*i])  
       { case Or: case ClosePar: case 0: /*or with next or end of pattern */
	     break;
	 default:  /* must be a concatenation */
	     e1 = malloc (sizeof(Tree));
	     e1->type = CONC;
	     e1->e1 = e; e = e1;
	     e->e2 = parseConc (pat,i,pos);
	     if (e->e2 == NULL) { freeTree(e); return NULL; }
	     e->eps = e->e1->eps && e->e2->eps;
	     e->maskPos = e->firstPos = e->lastPos = NULL;
	     break;
       }
    return e;
  }

static Tree *parseOr (byte *pat, int *i, Tree **pos)

	/* parses a disjunction */

  { Tree *e,*e1;
    e = parseConc (pat,i,pos); if (e == NULL) return NULL;
    switch (pat[*i])  
       { case ClosePar: case 0:  /* end of pattern */
	     break;
	 default:  /* must be an or */
	     e1 = malloc (sizeof(Tree));
	     e1->type = OOR;
	     e1->e1 = e; e = e1;
	     (*i)++;  /* skip the | */
	     e->e2 = parseOr (pat,i,pos);
	     if (e->e2 == NULL) { freeTree(e); return NULL; }
	     e->eps = e->e1->eps || e->e2->eps;
	     e->maskPos = e->firstPos = e->lastPos = NULL;
	     break;
       }
    return e;
  }

static void simpFree (Tree *e, Tree **pos, int m)

   { int i;
     if (e == NULL) return;
     simpFree (e->e1,pos,m); simpFree (e->e2,pos,m);
     for (i=0;i<m;i++)
	 if (pos[i] == e) pos[i] = NULL;
   }
     
static Tree *simplify (Tree *e, Tree **pos, int m, bool first, bool last)

   { int i;
     Tree *aux = e;
     if (e->eps && (first || last))  /* should be just an epsilon */
	{ simpFree (e,pos,m);
          freeTree(e->e1); e->e1 = NULL;
	  freeTree(e->e2); e->e2 = NULL;
	  e->type = STR; 
	  return e;
	}
     switch (e->type)
	{ case STR:   /* nothing to do */
	     break;
	  case OOR:   /* or of two classes is a class, C|eps -> C? */
	     e->e1 = simplify (e->e1,pos,m,first,false);
	     e->e2 = simplify (e->e2,pos,m,false,last);
	     if ((e->e1->type == STR) && (e->e2->type == STR) &&
		 (e->e1->eps == e->e2->eps))
		{ for (i=0;i<m;i++)
		      if ((pos[i] == e->e1) || (pos[i] == e->e2)) pos[i] = e;
		  e->type = STR; e->eps = e->e1->eps;
		  free (e->e1); e->e1 = NULL;
		  free (e->e2); e->e2 = NULL;
		}
	     else if ((e->e1->type == STR) && e->e1->eps &&
		      (e->e2->type == STR) && !e->e2->eps)
		  { e->e1->type = QUESTION;
		    e->e1->e1 = e->e2;
		    e = e->e1; free (aux);
		  }
	     else if ((e->e2->type == STR) && e->e2->eps &&
		      (e->e1->type == STR) && !e->e1->eps)
		  { e->e2->type = QUESTION;
		    e->e2->e1 = e->e1;
		    e = e->e2; free (aux);
		  }
	     else if ((e->e2->type == STR) && e->e2->eps &&
		      (e->e1->type == STR) && e->e1->eps)
		  { freeTree(e->e1); e = e->e2; free (aux);
		  }
	     break;
	  case CONC: /* conc with epsilon can be removed */
	     e->e1 = simplify (e->e1,pos,m,first,false);
	     e->e2 = simplify (e->e2,pos,m,false,last);
	     if ((e->e1->type == STR) && e->e1->eps)
		{ freeTree (e->e1); e = e->e2; free (aux); }
	     else if ((e->e2->type == STR) && e->e2->eps)
		{ freeTree (e->e2); e = e->e1; free (aux); }
	     break;
	  case STAR: /* ** = *, ?* = *, +* = *, eps* = eps */
	     e->e1 = simplify (e->e1,pos,m,first,last);
	     if ((e->e1->type == PLUS) || (e->e1->type == QUESTION) ||
		 (e->e1->type == STAR))
		{ e = e->e1; e->type = STAR; e->eps = true; free (aux); }
	     else if ((e->e1->type == STR) && e->e1->eps)
		{ e = e->e1; free (aux); }
	     break;
	  case QUESTION: /* *? = *, ?? = ?, +? = *, eps? = eps */
	     e->e1 = simplify (e->e1,pos,m,first,last);
	     if ((e->e1->type == STAR) || (e->e1->type == PLUS))
		{ e = e->e1; e->type = STAR; e->eps = true; free (aux); }
	     else if (e->e1->type == QUESTION)
		{ e = e->e1; free (aux); }
	     else if ((e->e1->type == STR) && e->e1->eps)
		{ e = e->e1; free (aux); }
	     break;
	  case PLUS: /* remove the + at beginning and ending
			    *+ = *, ?+ = *, ++ = +, eps+ = eps */
	     e->e1 = simplify (e->e1,pos,m,first,last);
	     if (first || last)
		{ e = e->e1; free (aux);
		}
	     else if ((e->e1->type == STAR) || (e->e1->type == QUESTION))
		{ e = e->e1; e->type = STAR; e->eps = true; free (aux); }
	     else if (e->e1->type == PLUS)
		{ e = e->e1; free (aux); }
	     else if ((e->e1->type == STR) && e->e1->eps)
		{ e = e->e1; free (aux); }
	     break;
	}
     return e;
   }

	/* determines the class of pattern for a subset active (all if NULL) */

static int detClass1 (Tree *e, Mask active, Mask tmp, int m)

   { int class1,class2;
     if (active && ISZERO(AND(COPY(tmp,active,m),e->maskPos,m),m)) 
	return SIMPLE;
     switch (e->type)
       { case STR:
	    return SIMPLE;
	 case QUESTION: case STAR: case PLUS:
	    class1 = detClass1 (e->e1,active,tmp,m);
	    if (class1 == SIMPLE) class1 = EXTENDED;
	    if (e->e1->type != STR) class1 = REGULAR;
	    return class1;
	 case OOR:
	    return REGULAR;
	 case CONC:
	    class1 = detClass1 (e->e1,active,tmp,m);
	    class2 = detClass1 (e->e2,active,tmp,m);
	    return max(class1,class2);
       }
   }

int detClass (Tree *e, Mask active, int m)

   { Mask tmp = createMask(m);
     int class = detClass1 (e,active,tmp,m);
     free (tmp);
     return class;
   }

static Tree *parseLiteral (byte *pat, int *i, Tree **pos)

	/* parses a literal pattern */

   { Tree *e,*r;
     e = malloc (sizeof(Tree));
     e->type = STR; e->eps = false;
     e->e1 = e->e2 = NULL;
     e->maskPos = e->firstPos = e->lastPos = NULL;
     pos[*i] = e; (*i)++;
     if (!pat[1]) r = e; /* one letter */
     else
	{ r = malloc (sizeof(Tree));
	  r->type = CONC; r->eps = false;
	  r->e1 = e; r->e2 = parseLiteral (pat+1,i,pos);
	}
     return r;
   }

void setMaskPos (Tree *e, int L)

        /* compute maskPos (positions of subexpression) */

  { switch (e->type)
       { case STR:
            e->maskPos = ZERO(createMask(L),L);
            if (!e->eps) SET(e->maskPos,e->pos);
            break;
         case STAR: case PLUS: case QUESTION:
            setMaskPos(e->e1,L);
            e->maskPos = COPY(createMask(L),e->e1->maskPos,L);
            break;
         case OOR: case CONC:
            setMaskPos(e->e1,L); setMaskPos(e->e2,L);
            e->maskPos =OR(COPY(createMask(L),e->e1->maskPos,L),
                           e->e2->maskPos,L);
            break;
       }
  }

Tree *parse (byte *pat, int m, Tree **pos)

	/* captures the structure of the regular expression and
	   maps each character to the place where it has to be put. It then
	   performs the optimizations and returns the real type of the search
	   pattern (simple, extended or regular) and the real number of bits
	   to store it. It returns NULL if there is a syntax error. */

  { Tree *e;
    int i;
    bool fl = !OptWholeWord && !OptWholeRecord;

    if (m==0) return NULL;  /* void pattern */
    for (i=0;i<m;i++) pos[i] = NULL;
    i = 0;
	/* parse the expression */
    if (OptLiteral) return parseLiteral (pat,&i,pos);
    e = parseOr (pat,&i,pos); if (e == NULL) return NULL;
	/* simplify the expression */
    e = simplify (e,pos,m,fl,fl);
    return e;
  }

