
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

	/* Search an extended pattern in a text */

#include "extended.h"
#include "simple.h"

#define BF 0.70  /* bwd has to be at most this to prefer it over fwd */

   	/* for reasons explained later, this must ensure that it never
   	   recommends a subpattern that begins or ends with ?s */

double extendedFindBest (Mask *B, Mask *S, Mask A, int m, int K,
		         int *wlen, int *beg, int *end)

	/* Finds a best factor to search from *beg to *end-1, and
	   puts in *wlen the window length (0 if it recommends 
	   forward scanning). K is a number of chars always assumed
	   to be scanned. Returns the avg number of chars inspected */

   { double *prob;    /* prob[i] = prob matching pattern position i */
     double *sprob;   /* sprob[i] = prob staying at pattern position i */
     double *pprob;   /* pprob[i,j,l] = prob of factor of length l 
			 starting at i and not surpassing character j */
     double *mprob;   /* mprob[i,j,l] = prob of any factor of length l
		         starting at >=i and not surpassing character j */
     int *lastl;      /* lastl[j] = up to which l is pprob[i,j,l]
			 computed (for fixed i) */
     int i,c,j,k,l,t,L;
     double cost,best;
     int m1,m2;

		/* first, load probabilities, O(m) time */

     prob = malloc (m * sizeof (double));
     sprob = malloc (m * sizeof (double));
     for (i=0; i<m; i++) 
       { prob[i] = sprob[i] = 0.0;
	 for (c=0; c<256; c++)
	     { if (ISSET(B[c],i)) prob[i] += letterProb[c];
	       if (ISSET(S[c],i)) sprob[i] += letterProb[c];
	     }
       }

     		/* second, prepare for lazy filling of pprob */

     pprob = malloc ((m+1) * m * (m+1) * sizeof (double));
     mprob = malloc ((m+1) * m * (m+1) * sizeof (double));
     m1= m+1; m2 = m * m1;
     lastl = malloc (m * sizeof (int));
     for (j=0;j<m;j++) 
	 { lastl[j] = 0;
	   for (i=0;i<=j;i++) pprob[i*m2+j*m1+0] = mprob[i*m2+j*m1+0] = 1.0;
	   pprob[(j+1)*m2+j*m1+0] = mprob[(j+1)*m2+j*m1+0] = 0.0;
	 }

		/* third, find best factor. this takes O(m^3) in the
		   worst case but should be close to O(m^2 log m) on
		   average */
    
     best = BF; *beg = *end = *wlen = 0;
     for (i=0; i<m; i++)
       { L = 0;
	 for (j=i+1; j<=m; j++)  /* analyze factor[i,j-1] */
	     { if (j-i > W) continue;
	       if (!ISSET(A,j-1)) L++;
               if (L <= 2*K) continue;
	       cost = K+1.0;
	       for (l=1; l<=L; l++)  /* avg length of scanned window */
		 { if ((cost>=L-K+1) || (cost/(L-K-cost+1) >= best))
		       break;  /* cannot win */
		   if (lastl[j-1] < l)  /* fill a row (always one) */
		      { int ii;
			pprob[j*m2+(j-1)*m1+l] = mprob[j*m2+(j-1)*m1+l] = 0.0;
			for (ii=j-1;ii>=0;ii--)
			   { pprob[ii*m2+(j-1)*m1+l] =
			        prob[ii] * pprob[(ii+1)*m2+(j-1)*m1+(l-1)] +
			        sprob[ii] * pprob[ii*m2+(j-1)*m1+(l-1)] +
			        (ISSET(A,ii) ? pprob[(ii+1)*m2+(j-1)*m1+l] : 0);
			     if (pprob[ii*m2+(j-1)*m1+l] > 1.0)
				pprob[ii*m2+(j-1)*m1+l] = 1.0;
			     mprob[ii*m2+(j-1)*m1+l] =
				           1-(1-pprob[ii*m2+(j-1)*m1+l])*
				             (1-mprob[(ii+1)*m2+(j-1)*m1+l]);
			   }
			lastl[j-1] = l;
		      }
		   cost += mprob[i*m2+(j-1)*m1+l];
		 }
	       if ((cost < L-K+1) && (cost/(L-K-cost+1) < best))
		  { best = cost/(L-K-cost+1);
		    *beg = i; *end = j; *wlen = L;
		  }
	     }
       }
     free (prob); free (sprob);
     free (pprob); free (mprob); free (lastl);

		/* final checks and possible switch to fwd scan */

     if (*wlen > 0)
	{       /* a check for safety: haven't I included ?s at
     		   the beginning or the end; right? */
          while ((*beg < *end) && ISSET(A,*beg)) (*beg)++;
          while ((*beg < *end) && ISSET(A,*end-1)) (*end)--;
	  if (*beg == *end) *wlen = 0;
	}
     if (*wlen == 0) /* fwd scan, this proc guarantees a legal subpattern */
	{ *end = m; if (*end > W) *end = W;
          while (ISSET(A,*end-1)) (*end)--;
	  best = 1.0;
	}
     return best;
   }

        /* receives L in zero and returns there the number of bits to hold
           the pattern */

void extendedLength (Tree *tree, int *L)

   { switch (tree->type)
        { case STR:
             tree->pos = (*L)++;
             break;
          case CONC:
             extendedLength(tree->e1,L);
             extendedLength(tree->e2,L);
             break;
	  case PLUS: case STAR: case QUESTION:
	     extendedLength(tree->e1,L);
             break;
        }
   }

static void extendedTreeLoad (Tree *tree, Mask *B, Mask *S, Mask A)

                /* reads masks S,A from pattern */

   { int c,l;
     switch (tree->type)
	{ case STR:
	    break;
	  case PLUS: /* add in S what we did to B */
	    l = tree->e1->pos;
	    for (c=0; c<256; c++) 
		if (ISSET(B[c],l)) SET(S[c],l);
	    break;
	  case QUESTION:  /* mark S */
	    SET(A,tree->e1->pos);
	    break;
	  case STAR:  /* treat it as +? */
	    l = tree->e1->pos;
	    for (c=0; c<256; c++) 
		if (ISSET(B[c],l)) SET(S[c],l);
	    SET(A,l);
	    break;
	  case CONC: /* recurse */
	    extendedTreeLoad (tree->e1,B,S,A);
	    extendedTreeLoad (tree->e2,B,S,A);
	    break;
	}
   }

void extendedLoadMasks (byte *pat, int m, int L, Tree *tree, Tree **pos,
			Mask *B, Mask *S, Mask *A)

		/* load masks B,S,A from pattern pat */

   { int i,c,l;
     B[0] = createMasks (256,L);
     S[0] = createMasks (256,L);
     for (c=0; c<256; c++)
	 { B[c] = ZERO(B[0]+c*maskSize(L),L);
	   S[c] = ZERO(S[0]+c*maskSize(L),L);
	 }
     (*A) = ZERO(createMask(L),L);

		/* first load B */
     i = 0;
     while (i<m)
        if (pos[i]) getAclass (pat,&i,B,pos[i]->pos);
        else i++;
     if (OptRecChar != -1) ZERO(B[OptRecChar],L); 

		/* now load S and A where appropriate */

     extendedTreeLoad (tree,B,S,*A);
  }

void extendedLoadVerif (int m, Mask *B, Mask *S, Mask A, int base, int sign,
		       Mask *BV, Mask *SV, Mask *IV, Mask *FV, Mask *AV,
		       Mask *DiV)

			/* loads the fwd or bwd masks */

   { int c,l;
     bool firstF;
     BV[0] = createMasks (256,m);
     SV[0] = createMasks (256,m);
     for (c=0; c<256; c++)
	 { BV[c] = m ? ZERO(BV[0]+c*maskSize(m),m) : NULL;
	   SV[c] = m ? ZERO(SV[0]+c*maskSize(m),m) : NULL;
	 }
     *IV = ZERO(createMask(m),m);
     *FV = ZERO(createMask(m),m);
     *AV = ZERO(createMask(m),m);
     *DiV = ZERO(createMask(m),m);

     firstF = false;
     for (l=0; l<m; l++)
         { for (c=0; c<256; c++) 
	     { if (ISSET(B[c],base+l*sign)) SET(BV[c],l);
	       if (ISSET(S[c],base+l*sign)) SET(SV[c],l);
	     }
	   if (ISSET(A,base+l*sign))
	      { if (l > 0)
		   if (!ISSET(*FV,l-1))
		      { SET(*IV,l-1); SET(*FV,l);
		        firstF = true;
		      }
	           else   /* not the first of a seq of ?'s */
		      { CLEAR(*FV,l-1); SET(*FV,l);
		      }
	        if (firstF) SET(*AV,l);
		else SET(*DiV,l);
	      }
	 }
   }

void extendedFreeScan (extendedScanData *scan)

   { free (scan);
   }

extendedScanData *extendedLoadFast (int wlen, Mask *B, Mask *S, Mask A, 
				    int beg, int end)

		/* loads masks for fast bwd/fwd */

   { int l,c;
     int base,sign,tbase;
     extendedScanData *scan = malloc (sizeof(extendedScanData));
     for (c=0; c<256; c++) scan->B[c] = scan->S[c] = ZEROS;
     scan->I = scan->F = scan->A = ZEROS;
     scan->m = end-beg;
     scan->wlen = wlen;
     if (wlen)
	  { base = end-1; sign = -1; tbase = W-scan->m; }
     else { base = beg; sign = +1; tbase = 0; }
     for (l=0; l<scan->m; l++)
	{ for (c=0; c<256; c++) 
	     { if (ISSET(B[c],base+sign*l)) scan->B[c] |= ONE << (tbase+l); 
	       if (ISSET(S[c],base+sign*l)) scan->S[c] |= ONE << (tbase+l);
	     }
	  if (ISSET(A,base+sign*l))
	     { scan->A |= ONE << (tbase+l);
	       if (!(scan->F & (ONE << (tbase+l-1))))
		  { scan->I |= ONE << (tbase+l-1);
		    scan->F |= ONE << (tbase+l);
		  }
	       else   /* not the first of a seq of ?'s */
		  { scan->F &= ~ (ONE << (tbase+l-1));
		    scan->F |= ONE << (tbase+l);
		  }
	     }
       }
     
     return scan;
   }

   	/* the real parser will trim extra ?s and *s at the extremes of pat */

extendedData *extendedPreproc (byte *pat, Tree *tree, Tree **pos)

	/* Preprocesses pat and creates an extendedData structure
	   for searching */

   { extendedData *P = malloc (sizeof(extendedData));
     Mask B[256],S[256],A,active;
     int i,l,c,L,wlen,beg,end;

	/* first, compute the length of the pattern and
	   check that it has extended syntax */

     L = 0;
     extendedLength (tree,&L);
     setMaskPos (tree,L);


	/* second, repeat the process creating bit masks */

     extendedLoadMasks (pat,strlen(pat),L,tree,pos,B,S,&A);

	/* third, find the best factor to search the pattern */

     extendedFindBest (B,S,A,L,0,&wlen,&beg,&end);
     P->type = (wlen == 0) ? FWD: BWD;

	/* fourth, create the data structure according to the above
	   decision */

     P->mbwdV = wlen ? beg : end;
     extendedLoadVerif (P->mbwdV,B,S,A,P->mbwdV-1,-1, P->BbwdV,P->SbwdV,
		        &P->IbwdV,&P->FbwdV,&P->AbwdV,&P->DibwdV);
     P->mfwdV = wlen ? L-beg : L-end;
     extendedLoadVerif (P->mfwdV,B,S,A,L-P->mfwdV,+1, P->BfwdV,P->SfwdV,
		     	&P->IfwdV,&P->FfwdV,&P->AfwdV,&P->DifwdV);

     	/* when selecting a subpattern it is possible that an I bit
     	   dos not have its F in the subpattern (yes, this means a bad
	   optimization algorithm). but even with the missing F bit the
	   I bit will work ok with the subtraction. an F bit without its
	   I bit will not work properly except if we do some extra work.
	   since this is a bad idea we ensure that the optimization will
	   never recommend such a pattern. finally, an I bit at the last 
	   position is meaningless and can be discarded (and it is not an
	   error in the optimization scheme) */
	  
     active = ZERO(createMask(L),L);
     for (i=beg;i<end;i++) SET(active,i);
     if (detClass(tree,active,L) == EXTENDED)
        { P->scanData = extendedLoadFast (wlen,B,S,A,beg,end);
	  P->scanText = (scanProc)extendedScan;
	  P->scanFree = (freeScanProc)extendedFreeScan;
	}
     else 
	{ P->scanData = simpleLoadFast (wlen?true:false,B,beg,end);
	  P->scanText = (scanProc)simpleScan;
	  P->scanFree = (freeScanProc)simpleFreeScan;
	}

	/* finally, free temporary structures and return P */

     free (B[0]); free (S[0]); free (A);
     free (active);
     P->DV = createMask (L);
     return P;
   }

void extendedFree (extendedData *P)

	/* Frees P */

   { P->scanFree (P->scanData);
     free (P->BbwdV[0]); free (P->SbwdV[0]);
     free (P->IbwdV); free (P->FbwdV); free (P->AbwdV);
     free (P->DibwdV);
     free (P->BfwdV[0]); free (P->SfwdV[0]);
     free (P->IfwdV); free (P->FfwdV); free (P->AfwdV);
     free (P->DifwdV);
     free (P->DV);
     free (P);
   }


static byte *dirCheck (byte *ptr, byte *top, int inc, 
		       int m, Mask *B, Mask *S,
		       Mask I, Mask F, Mask A,
		       Mask Di, Mask D)

		/* directional check, forward or backward.
		   returns NULL if not found, last pos read if found */
		/* assumes that ptr needs to be shifted before reading */

		/* note that this subpattern CAN have ?s or *s at the
		   beginning (but still not at the end) */

   { mask Df;			/* masks */ 
     mask X,E;			/* 1 then 0 */
     Mask nB,nS;
     int r,pos;
     mask carry,ohigh,nhigh;
     bool alive;
     int n;
     
        /* start by removing case m = 0 */

     if (m == 0)
	{ if ((inc == -1) && recCheckLeftContext(ptr,top)) return ptr;
	  if ((inc == +1) && recCheckRightContext(ptr+1,top+1)) return ptr;
	  return NULL;
	} 

     	/* we first account for initial ?s. if the first F has no
     	   matching I we put 1s up to its position - 1 */

     n = (m+W-1)/W;
     COPY (D,Di,m);
     X = ONE;
     E = ONE << ((m-1)%W);

		/* multimasks not used for efficiency */

     while (true)
	{ if (E & D[n-1])  /* match, check context */
	     { if ((inc == -1) && recCheckLeftContext(ptr,top)) return ptr;
	       if ((inc == +1) && recCheckRightContext(ptr+1,top+1)) return ptr;
	     }
	  if (ptr == top) return NULL;  /* end of buffer */
	  ptr += inc;
	  nB = B[*ptr]; nS = S[*ptr];
	     /* D = (((D << 1) | X) & B[c]) | (D & S[c]); */
	  ohigh = X; X = ZEROS;
          alive = false;
          for (r=0; r<n; r++) 
              { nhigh = (D[r] >> (W-1)) & ONE;
	        D[r] = (((D[r] << 1) | ohigh) & nB[r]) | (D[r] & nS[r]);
	        ohigh = nhigh;
		if (D[r]) alive = true;
	      }
	  if (!alive) return NULL;  /* automaton died */
             /* Df = D | F; D |= A & ((~(Df-I)) ^ Df); */
	  carry = ZEROS;
          for (r=0; r<n; r++) 
              { Df = D[r] | F[r];
		D[r] |= A[r] & ((~(Df-carry-I[r])) ^ Df);
		carry = (Df < carry + I[r]) || (carry + I[r] < carry);
	      }
        }
   }

static bool checkMatch (extendedData *P, byte *pos, byte **beg, byte **end)

	/* Checks that an initially promising match starting at pos
	   (if P->bwd) or ending at pos-1 (if !P->bwd)
	   corresponds or not to a complete match. The buffer limits
	   are *beg and *end, and in case of successful match we return
	   there the boundaries of the match found. */

   { byte *rbeg,*rend,*obeg,*oend;

        /* start by knowing my surrounding record */

     if (P->type == FWD)  /* pos-1 is last position read */
        { recGetRecord (pos-1,*beg,*end,&rbeg,&rend,&obeg,&oend);
          if ((rbeg > pos-1) || (rend <= pos-1)) return false; /* overlaps */
        }
     else
        { recGetRecord (pos,*beg,*end,&rbeg,&rend,&obeg,&oend);
          if ((rbeg > pos) || (rend <= pos)) return false; /* overlaps */
        }   

	/* first part: backward check */

     if (dirCheck (pos,rbeg,-1,P->mbwdV,P->BbwdV,P->SbwdV,
		   P->IbwdV,P->FbwdV,P->AbwdV,P->DibwdV,P->DV) == NULL)
	return false;

	/* second part: forward check */

     if (dirCheck (pos-1,rend-1,+1,P->mfwdV,P->BfwdV,P->SfwdV,
		   P->IfwdV,P->FfwdV,P->AfwdV,P->DifwdV,P->DV) == NULL)
	return false;

	/* match found */

     *beg = obeg;
     *end = oend;
     return true;
   }

bool extendedScan (byte **beg, byte **end, checkProc checkMatch, void *P,
		   extendedScanData *scan)
 
        /* Scans from *beg to *end-1 for P. Returns if it could find
           it. In case of returning true, *beg and *end are set to limit
           the first occurrence of P. It requires
             - a procedure
                bool checkMatch (P,byte *pos, byte **beg, byte **end)
               which whether there is a match around pos and returns in
               *beg and *end the limits. "Around" means that pos is the initial               position of the search pattern given to extendedScan.
             - tables of masks B,S,I,F,A for the pattern
             - length M of the search pattern
	     - window length m (zero => fwd scan)
         */                                                                    

   { register byte *pos,*top;	/* current and final text pointers */
     register mask D,Df;	/* D masks for scanning */
     register mask E;
     register int j;		/* window index */
     register int c;
     register mask *B = scan->B;
     register mask *S = scan->S;
     register mask I = scan->I;
     register mask F = scan->F;
     register mask A = scan->A;
     register int m = scan->wlen;

     pos = *beg; top = *end;
     
     if (m > 0)    /* backward scanning */
        { pos--;
	  top -= m;
	  while (pos < top)
             { D = B[pos[m]]; if (!D) { pos += m; continue; } /* skip loop */
	       j = m;
	       do { if (--j == 0) 
	               {  /* a factor of the pattern has matched, 
			     if it is a prefix, check all */
	                 if ((D & (ONE<<(W-1))) && checkMatch (P,pos+1,beg,end))
			    return true;
	 		  /* else shift in 1 */
	                 j = 1;
			 break;
	               }
	            Df = D | F;	/* expand questionmarks */
		    D |= A & ((~(Df-I)) ^ Df);
                    c = pos[j]; /* take new character */
	            D = ((D << 1) & B[c]) | (D & S[c]);  /* stars */
		  }
	       while (D);
	       pos += j;
             }
        }
     else	/* forward scanning */
        { E = ONE << (scan->m-1);
          while (pos < top)
            { D = ZEROS;
	      while (pos < top)
	        { c = *pos++;
                  D = (((D << 1) | ONE) & B[c]) | (D & S[c]);  /* stars */
	          Df = D | F; /* expand questionmarks */
	          D |= A & ((~(Df-I)) ^ Df);
                  if (D & E)
	             {  /* selected part of pattern has matched, check all */
	               if (checkMatch (P,pos,beg,end)) return true;
	             }
                }
            }
        }
     return false;
   }

bool extendedSearch (byte **beg, byte **end, extendedData *P)

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

	/* it does not work if there are ?s at the beginning or the end,
	   so this has to be excluded by the optimizer */

   { return P->scanText (beg,end,(checkProc)checkMatch,P,P->scanData);
   }
 
