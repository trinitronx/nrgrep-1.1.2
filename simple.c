
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

	/* Search a simple pattern in a text */

#include "simple.h"

#define BF 0.80  /* bwd has to be at most this to prefer it over fwd */

double simpleFindBest (Mask *B, int m, int K,
		       bool *bwd, int *beg, int *end)

	/* Finds a best factor to search from *beg to *end-1, and
	   recommends in *bwd a bwd o fwd search of it. Assumes that
	   K extra characters are read in any window. Returns avg
	   work per text character. */

   { double *prob;   /* prob[i] = prob matching pattern pos i */
     double *pprob;  /* pprob[i,l] = prob matching P[i..i+l-1] */
     double *mprob;  /* for fixed i, mprob[l-1] = prob of any subpattern of
			length l in [i..x] */
     int *lprob;     /* value of x */
     int i,c,j,k,l,t;
     double best,cost;

		/* first, load probabilities */

     prob = malloc (m * sizeof (double));
     for (i=0; i<m; i++) 
       { prob[i] = 0.0;
	 for (c=0; c<256; c++)
	     if (ISSET(B[c],i)) prob[i] += letterProb[c];
       }

     		/* second, compute matching prob of a factor of length l
     		   starting at position k. This takes O(m^2) time. */

     pprob = malloc ((m+1) * (m+1) * sizeof (double));
     pprob[(m+1)*m+0] = 1.0;
     for (l=1;l<=m;l++) pprob[(m+1)*m+l] = 0.0;
     for (i=m-1;i>=0;i--)
         { pprob[(m+1)*i+0] = 1.0;
	   for (l=1;l<=m;l++)
               pprob[(m+1)*i+l] = prob[i] * pprob[(m+1)*(i+1)+(l-1)];
	 }
     free (prob);
	
		/* third, find best factor. this takes O(m^3) in the
		   worst case but is closer to O(m^2 log m) on average */
    
     best = BF; *beg = *end = 0;
     mprob = malloc (m * sizeof (double));
     lprob = malloc (m * sizeof (int));
     for (i=0; i<m; i++)
       {    /* initialize vector of max probs of length l from i to j-l */
	 for (l=1; l<=m; l++)  /* avg length of scanned window */
	     { mprob[l-1] = 0.0;/* max prob length l in [i..something] */
	       lprob[l-1] = i+l-2;  /* value of "something" */
	     }
	    /* now find the best combination */
         for (j=i+K+1; (j<=m) && (j-i<=W); j++)  /* analyze factor[i,j-1] */
	     { cost = K+1;
	       for (l=1; l<=j-i; l++)  /* avg length of scanned window */
		 { if ((cost>=j-i-K+1) || (cost/(j-i-K-cost+1) >= best))
		      break; /* cannot win */
		   for (k=lprob[l-1]+1; k<=j; k++)
		       mprob[l-1] = 1-(1-mprob[l-1])*(1-pprob[(m+1)*(k-l+1)+l]);
		   lprob[l-1] = j;
		   cost += mprob[l-1];
		 }
	       if ((cost < j-i-K+1) && (cost/(j-i-K-cost+1) < best))
		  { best = cost/(j-i-K-cost+1);
		    *beg = i; *end = j;
		  }
	     }
       }
     free (pprob);
     free (mprob);
     free (lprob);

     		/* just a consistency check: didn't we set end-beg <= k+1 ? */

     if (*end-*beg <= K+1) { *beg = *end = 0; }

		/* finally, determine fwd or bwd */

     *bwd = (*end != 0);  /* fwd seems to be the best choice */
     if (!*bwd)
	{ *end = m; if (*end > W) *end = W;
	}

     return best < BF ? best : 1.0;
   }

	/* receives L in zero and returns there the number of bits to hold
	   the pattern */

void simpleLength (Tree *tree, int *L)

   { switch (tree->type)
	{ case STR: 
	     tree->pos = (*L)++;
	     break;
	  case CONC: 
	     simpleLength(tree->e1,L);
	     simpleLength(tree->e2,L);
	     break;
	}
   }

void simpleLoadMasks (byte *pat, int m, int L, Mask *B, Tree **pos)

		/* reads mask B from pattern */

   { int i,c,l;
     B[0] = createMasks (256,L);
     for (c=0; c<256; c++)
	 B[c] = ZERO(B[0]+c*maskSize(L),L);

     i = 0;
     while (i<m)
	if (pos[i])
	   { if (OptLiteral) { SET(B[pat[i]],pos[i]->pos); i++; }
	     else getAclass (pat,&i,B,pos[i]->pos);
	   }
	else i++;
     if (OptRecChar != -1) ZERO(B[OptRecChar],L);
   }

void simpleLoadVerif (int m, Mask *B, int base, int sign, Mask *BV)

		/* loads bwd or fwd verification masks */

   { int c,l;
     BV[0] = createMasks (256,m);
     for (c=0; c<256; c++)
	 BV[c] = m ? ZERO(BV[0]+c*maskSize(m),m) : NULL;
     for (l=0; l<m; l++)
	 for (c=0; c<256; c++) 
	     if (ISSET(B[c],base+l*sign)) SET(BV[c],l);
   }

void simpleFreeScan (simpleScanData *scan)

   { free (scan);
   }

                /* loads masks for fast bwd/fwd */

simpleScanData *simpleLoadFast (bool bwd, Mask *B, int beg, int end)

   { int c,l;
     simpleScanData *scan = malloc (sizeof(simpleScanData));
     scan->m = end-beg;
     if (bwd) 
	{ scan->bwd = true;
          for (c=0; c<256; c++) scan->B[c] = ZEROS;
          for (l=0; l<scan->m; l++)
	     for (c=0; c<256; c++) 
	        if (ISSET(B[c],end-1-l)) scan->B[c] |= ONE << (W-scan->m+l); 
	}
     else
	{ scan->bwd = false;
          for (c=0; c<256; c++) 
	     scan->B[c] = (scan->m == W) ? ONES : ((ONE << scan->m)-1);
          for (l=0; l<scan->m; l++)
	     for (c=0; c<256; c++) 
	        if (ISSET(B[c],beg+l)) scan->B[c] &= ~ (ONE << l); 
	}
     return scan;
   }

simpleData *simplePreproc (byte *pat, Tree *tree, Tree **pos)

	/* Preprocesses pat and creates a simpleData structure
	   for searching */

   { simpleData *P = malloc (sizeof(simpleData));
     Mask B[256];
     int c,L;
     bool bwd;
     int beg,end;

	/* first, compute the length of the pattern in bits */

     L = 0;
     simpleLength (tree,&L);

	/* second, create bit masks */

     simpleLoadMasks (pat,strlen(pat),L,B,pos);

	/* third, find the best factor to search the pattern */

     simpleFindBest (B,L,0,&bwd,&beg,&end);

	/* fourth, create the data structure according to the above
	   decision */

     P->m = end-beg;
     P->mbwdV = beg;
     simpleLoadVerif (P->mbwdV,B,beg-1,-1,P->BbwdV);
     P->mfwdV = L-end;
     simpleLoadVerif (P->mfwdV,B,end,+1,P->BfwdV);

     P->scanData = simpleLoadFast (bwd,B,beg,end);
     P->scanText = (scanProc)simpleScan;
     P->scanFree = (freeScanProc)simpleFreeScan;

	/* finally, free temporary structures and return P */

     free (B[0]);
     P->raw = false;
     return P;
   }

void simpleFree (simpleData *P)

	/* Frees P */

   { P->scanFree (P->scanData);
     free (P->BbwdV[0]);
     free (P->BfwdV[0]);
     free (P);
   }

static bool checkMatch (simpleData *P, byte *pos, byte **beg, byte **end)

	/* Checks that an initially promising match starting at pos
	   corresponds or not to a complete match. The buffer limits
	   are *beg and *end, and in case of successful match we return
	   there the boundaries of the match found. */

   { Mask *Bf = P->BfwdV;	/* B table for fwd verif */
     Mask *Bb = P->BbwdV;	/* B table for bwd verif */
     byte *ptrb,*ptrf;		/* ptr to match of selected subpattern */
     int mb,mf;			/* verif length */
     int jb,jf;			/* window position */
     byte *text,*top;		/* buffer limits */
     byte *obeg,*oend;		/* record to return */

     if (P->raw)  /* no records, no contexts */
        { text = *beg; top = *end; }
     else 
	{ recGetRecord (pos,*beg,*end,&text,&top,&obeg,&oend);
	  if ((pos - P->mbwdV < text) || (pos + P->m + P->mfwdV > top))
	     return false;  /* overlaps with record limits */
	}

	/* first part: backward check */

     jb = 0;
     ptrb = pos - 1;
     mb = P->mbwdV;
     while ((jb < mb) && (ptrb-jb >= text))
	{ if (!ISSET(Bb[ptrb[-jb]],jb)) break;
	  jb++;
	}
     if (jb < mb) return false;  /* no match */

	/* second part: forward check */

     jf = 0;
     ptrf = pos + P->m;
     mf = P->mfwdV;
     while ((jf < mf) && (ptrf+jf < top))
	{ if (!ISSET(Bf[ptrf[jf]],jf)) break;
	  jf++;
	}
     if (jf < mf) return false;  /* no match */

        /* third part: check context */

     if (!P->raw)
	if (!recCheckLeftContext (ptrb-jb+1,text) ||
	    !recCheckRightContext (ptrf+jf,top))
	return false;  /* context failed */

	/* match found */

     if (P->raw)
          { *beg = ptrb-jb+1; *end = ptrf+jf; }
     else { *beg = obeg; *end = oend; }
     return true;
   }

bool simpleScan (byte **beg, byte **end, checkProc checkMatch, void *P,
		 simpleScanData *scan)

	/* Scans from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P. It requires
	     - a procedure 
	        bool checkMatch (P,byte *pos, byte **beg, byte **end)
	       which whether there is a match around pos and returns in
               *beg and *end the limits. "Around" means that pos is the initial
	       position of the search pattern given to simpleScan.
	     - a table of masks B for the pattern, negated if !bwd
	     - length m of the search pattern
	     - bwd, recommending bwd scanning (or otherwise fwd scanning)
	 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask D;		/* D mask for scanning */
     register int j;		/* window index */
     register mask *B = scan->B;
     register int m = scan->m;
     register mask E = ONE<<(m-1);/* E mask */

     pos = *beg; top = *end;
     
     if (scan->bwd)    /* backward scanning */
        { pos--;
	  top -= m;
	  while (pos < top)
             { D = B[pos[m]]; if (!D) { pos += m; continue; }  /* skip loop */
	       j = m;
               do D = (D << 1) & B[pos[--j]]; while (D);
               if (j == 0) 
	          {  /* selected part of pattern has matched, check the rest */
	            if (checkMatch (P,pos+1,beg,end)) return true;
	 		/* else shift in 1 */
	            j = 1;
	          }
	       pos += j;
             }
        }
     else	/* forward scanning */
        { D = ONES;
          while (pos < top)
            { D = (D << 1) | B[*pos++];
              if (!(D & E))
	         {  /* selected part of pattern has matched, check the rest */
	           if (checkMatch (P,pos-m,beg,end)) return true;
	         }
            }
        }
     return false;
   }

bool simpleSearch (byte **beg, byte **end, simpleData *P)

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

   { return P->scanText (beg,end,(checkProc)checkMatch,P,P->scanData);
   }

bool simpleRevSearch (byte **beg, byte **end, simpleData *P)

	/* Searches from *end-1 to *beg for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the last occurrence of P (first found in bwd scan) */
	/* Knows for sure that the pattern is simple */

   { simpleScanData *scan = P->scanData;
     register byte *pos,*top;	/* current and final text pointers */
     register mask D;		/* D mask for scanning */
     register mask E = ONE;	/* E mask */
     register mask *B = scan->B;/* B table for scanning */
     register int j;		/* window index */
     register int m = scan->m;	/* len of sel. pattern / window length */
     register mask I;		/* Initial mask */

     pos = *end-1; top = *beg;
     I = (m == W) ? ONES : ((ONE << m)-1);
     
     if (scan->bwd)    /* backward scanning */
        { pos++;
	  m = -m;
	  top -= m;
	  while (pos >= top)
             { D = B[pos[m]]; if (!D) { pos += m; continue; } /* skip loop */
               j = m;
               do D = (D >> 1) & B[pos[++j]]; while (D);
               if (j == 0) 
	          {  /* selected part of pattern has matched, check the rest */
	            if (checkMatch (P,pos+m,beg,end)) return true;
	 		/* else shift in 1 */
		    j = -1;
	          }
	       pos += j;
             }
        }
     else	/* forward scanning */
        { D = I;
          while (pos >= top)
            { D = (D >> 1) | B[*pos--];
              if (!(D & E))
	         {  /* selected part of pattern has matched, check the rest */
	           if (checkMatch (P,pos+1,beg,end)) return true;
	         }
            }
        }
     return false;
   }

