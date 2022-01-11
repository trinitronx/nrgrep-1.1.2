
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

	/* Search an extended pattern in a text allowing errors */

#include "eextended.h"
#include "esimple.h"
#include "parser.h"
#include "options.h"

#define BK 0.95  /* k+1 has to be at most this to prefer it over fwd/bwd */

static double findBestMulti (Mask *B, Mask *S, Mask A, int m,
	       		     int K, int *wlen, int *beg, int *end)

	/* Finds a best set of K+1 factors to search, puts their length
	   in *wlen and their initial/final positions in beg[i]/end[i]. 
	   If no good set exists it returns wlen=0
	*/

   { int tr = OptTransp ? 1 : 0;
     int L = min(m-tr*K,W)/(K+1);
     double *prob;  /* prob[i] = prob pattern letter i */
     double *sprob; /* sprob[i] = prob staying at pattern position i */
     double *eprob; /* eprob[i,j,l] = prob of factor of length l
		       starting at i and not surpassing character j */
     double *mprob; /* mprob[i,j,l] = prob of any factor of length l
		       starting at >=i and not surpassing character j */
     int *lastl;    /* lastl[j] = up to which l is eprob[i,j,l]
		       computed (for fixed i) */                             
     double *pcost; /* pcost[i,j] = average cost per window of factor [i,j-1] */
     double *mbest; /* mbest[i,k] = cost using best selection of k strings
		       starting at i. the length of the strings is fixed */
     int *ibest;    /* where to start the first strings to obtain the
		       mbest cost */
     int *reach;    /* reach[i,l] = last symbol in P reachable in l iters */
     int wl,i,c,j,k,l,t;
     double this,best,cost;
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

     		/* second, compute reach in O(m^2/k) time */

     reach = malloc ((L+1+tr) * m * sizeof(int));
     for (i=0;i<m;i++)
	 { for (l=0;l<=L+tr;l++)
	       { t = (l == 0) ? i : reach[(l-1)*m+i]+1;
		 if (t > m) t = m;
		 if (t > 0)
		    while ((t<m) && ISSET(A,t-1)) t++;
		 reach[l*m+i] = t;
	       }
         }

                /* third, prepare for lazy filling of eprob */
      
     eprob = malloc ((m+1) * m * (m+1) * sizeof (double));
     mprob = malloc ((m+1) * m * (m+1) * sizeof (double));
     m1= m+1; m2 = m * m1;
     lastl = malloc (m * sizeof (int));
     for (j=0;j<m;j++)
         { lastl[j] = 0;
           for (i=0;i<=j+1;i++) eprob[i*m2+j*m1+0] = mprob[i*m2+j*m1+0] = 1.0;
	 }
		     
     		/* fourth, compute pcost in O(m^3) time */

     pcost = malloc (m * m * sizeof (double));
     for (i=0; i<m; i++)
	{ for (j=i+1; j<=m; j++)
	     { if (j-i > W/(K+1)) { pcost[i*m+j] = 1.0; continue; }
		     /* above is simple heuristic to avoid exceeding W bits */
	       cost = 1.0;
               for (l=1; l<=j-i; l++)  /* length of scanned window */
	           { if (lastl[j-1] < l)  /* fill a row (always one) */
		        { int ii;
		          eprob[j*m2+(j-1)*m1+l] = mprob[j*m2+(j-1)*m1+l] =0.0;
		          for (ii=j-1;ii>=0;ii--)
			      { eprob[ii*m2+(j-1)*m1+l] =
		                   prob[ii] * eprob[(ii+1)*m2+(j-1)*m1+(l-1)] +
		                   sprob[ii] * eprob[ii*m2+(j-1)*m1+(l-1)] +
		                   (ISSET(A,ii)?eprob[(ii+1)*m2+(j-1)*m1+l]:0);
				if (eprob[ii*m2+(j-1)*m1+l] > 1.0)
				   eprob[ii*m2+(j-1)*m1+l] = 1.0;
				mprob[ii*m2+(j-1)*m1+l] = 1-
				         (1-eprob[ii*m2+(j-1)*m1+l])*
					 (1-mprob[(ii+1)*m2+(j-1)*m1+l]);
			      }
		          lastl[j-1] = l;
		        }
		     cost += mprob[i*m2+(j-1)*m1+l];
	           }
	       pcost[i*L+(j-1)] = cost;
	     }
	 }
     free (prob); free (sprob);
     free (eprob); free (mprob); free (lastl);

		/* fifth, find best factors. this takes O(m^2) time in
		   the worst case but the average should be closer to O(km) */
    
     mbest = malloc ((m+1) * (K+2) * sizeof (double));
     ibest = malloc ((m+1) * (K+2) * sizeof (int));
     *wlen = 0; best = BK;
     for (wl = L; wl > 1; wl--)
	{ if (best < 1/(double)wl) break;  /* cannot win */
	  for (i=0;i<=m;i++) mbest[i*(K+2)+(0)] = 0.0;
	  for (k=1;k<=K+1;k++) mbest[m*(K+2)+k] = 1.0;
	  for (k=1; k<=K+1; k++)
	     for (i=m-1; i>=0; i--)
                 { j = reach[wl*m+i];
		   cost = (pcost[i*L+(j-1)] < j-i+1) ?
			  pcost[i*L+(j-1)]/(j-i-pcost[i*L+(j-1)]+1) : 1.0;
		   if (cost > 1.0) cost = 1.0;
		   this = 1-(1-cost)*(1-mbest[reach[(wl+tr)*m+i]*(K+2)+(k-1)]);
		   ibest[i*(K+2)+k] = i;
	           if (mbest[(i+1)*(K+2)+k] < this)
		      { this = mbest[(i+1)*(K+2)+k];
			ibest[i*(K+2)+k] = ibest[(i+1)*(K+2)+k];
		      }
                   mbest[i*(K+2)+k] = this;
		 }
	  if (mbest[0*(K+2)+(K+1)] < best)  /* put the results in output */
	     { *wlen = wl;
	       i = 0;
	       for (k=K+1;k>=1;k--)
		  { beg[K+1-k] = ibest[i*(K+2)+k];
		    end[K+1-k] = reach[wl*m+beg[K+1-k]];
		    i = reach[(wl+tr)*m+beg[K+1-k]];
		  }
	       best = mbest[0*(K+2)+(K+1)];
	     }
	}
     free (reach); free (pcost);
     free (mbest); free (ibest);

     if (best >= BK) { *wlen = 0; return 1.0; }

     		/* safety check: didn't we select subpatterns starting or
     		   ending with A bits on, right? */

     for (i=0;i<=K;i++)
	 { while ((beg[i] < end[i]) && ISSET(A,beg[i])) beg[i]++;  
					/* note that this cannot alter wlen! */
	   while ((beg[i] < end[i]) && ISSET(A,end[i]-1)) end[i]--;
	   if (beg[i] == end[i])  /* something horrible has happened, give up */
	      { *wlen = 0; return 1.0; }
	 }

     return best;
   }

void eextendedFreeScan (eextendedScanData *scan)

   { extendedFreeScan(scan->edata);
     free (scan->V3);
     free (scan->ends);
     free (scan);
   }

eextendedScanData *eextendedLoadFast (int wlen, int k, int type, 
				      Mask *B, Mask *S, Mask A, 
			 	      int *beg, int *end)

   { int i,c,l,acc;
     eextendedScanData *scan = malloc (sizeof(eextendedScanData));
     scan->k = k;
     scan->type = type;
     if (type != KPLUS1)
	{ scan->edata = extendedLoadFast (wlen,B,S,A,*beg,*end);
	  scan->m = *end-*beg;
	  scan->V3 = malloc ((2*(k+1))*sizeof(mask));
	  scan->ends = NULL;
	}
     else   /* partition in k+1 */
        { scan->m = wlen;
	  scan->edata = malloc (sizeof(extendedScanData));
          scan->edata->wlen = wlen;
          scan->ends = malloc ((k+1) * sizeof(mask));
          scan->V3 = NULL;
          acc = 0;
          for (c=0; c<256; c++) scan->edata->B[c] = scan->edata->S[c] = ZEROS;
          scan->edata->I = scan->edata->F = scan->edata->A = ZEROS;
          for (i=0;i<=k;i++)
              { for (l=0; l<end[i]-beg[i]; l++)
                   { for (c=0; c<256; c++)
                        { if (ISSET(B[c],end[i]-1-l))
                            { scan->edata->B[c] |= ONE << (acc+l);
                              if (l > 0) scan->B0[c] |= ONE << (acc+l);
                            }
                          if (ISSET(S[c],end[i]-1-l)) scan->edata->S[c] |= ONE << (acc+l);
                        }
                     if (ISSET(A,end[i]-1-l))
                        { scan->edata->A |= ONE << (acc+l);
                          if (!(scan->edata->F & (ONE << (acc+l-1))))
                             { scan->edata->I |= ONE << (acc+l-1);
                               scan->edata->F |= ONE << (acc+l);
                             }
                          else   /* not the first of a seq of ?'s */
                             { scan->edata->F &= ~ (ONE << (acc+l-1));
                               scan->edata->F |= ONE << (acc+l);
                             }
                        }
                   }
                acc += end[i]-beg[i];
                scan->ends[i] = ONE << (acc-1);
              }
        }
     return scan;
   }

eextendedData *eextendedPreproc ( byte *pat, Tree *tree, Tree **pos, int k)

	/* Preprocesses pat and creates an eextendedData structure
	   for searching */

   { eextendedData *P = malloc (sizeof(eextendedData));
     Mask B[256],S[256],A,active;
     int i,c,c1,c2,l,L,wlen,wlen1;
     int *beg,*end,beg1,end1;
     double best1,best2;
     int K;

        /* first, compute the length of the pattern and
           check that it has extended syntax */

     L = 0;
     extendedLength (tree,&L);
     setMaskPos (tree,L);

        /* second, repeat the process creating bit masks */

     extendedLoadMasks (pat,strlen(pat),L,tree,pos,B,S,&A);

	/* third, find the best way to search the pattern */

     best1 = (k+1) * extendedFindBest (B,S,A,L,k,&wlen1,&beg1,&end1); 
     beg = malloc ((k+1)*sizeof(int));
     end = malloc ((k+1)*sizeof(int));
     best2 = findBestMulti (B,S,A,L,k,&wlen,beg,end); 
     if ((best1 <= best2) || !wlen)  /* we prefer suffix */
	{ K = 0; beg[0] = wlen1 ? beg1 : end1; end[0] = L-beg[0]; }
     else K = k;
	
	/* fourth, create the verification data, common proc */

     P->k = k;
     P->mbwdV = malloc ((K+1)*sizeof(int));
     P->BbwdV = malloc ((K+1)*sizeof(Mask*));
     P->SbwdV = malloc ((K+1)*sizeof(Mask*));
     P->FbwdV = malloc ((K+1)*sizeof(Mask));
     P->IbwdV = malloc ((K+1)*sizeof(Mask));
     P->AbwdV = malloc ((K+1)*sizeof(Mask));
     P->DibwdV = malloc ((K+1)*sizeof(Mask));
     P->mfwdV = malloc ((K+1)*sizeof(int));
     P->BfwdV = malloc ((K+1)*sizeof(Mask*));
     P->SfwdV = malloc ((K+1)*sizeof(Mask*));
     P->FfwdV = malloc ((K+1)*sizeof(Mask));
     P->IfwdV = malloc ((K+1)*sizeof(Mask));
     P->AfwdV = malloc ((K+1)*sizeof(Mask));
     P->DifwdV = malloc ((K+1)*sizeof(Mask));

     for (i=0;i<=K;i++) 
	 { P->mbwdV[i] = beg[i];
	   P->BbwdV[i] = malloc (256 * sizeof(Mask));
	   P->SbwdV[i] = malloc (256 * sizeof(Mask));
           extendedLoadVerif (P->mbwdV[i],B,S,A,P->mbwdV[i]-1,-1,
	                      P->BbwdV[i],P->SbwdV[i],&P->IbwdV[i],
			      &P->FbwdV[i],&P->AbwdV[i],&P->DibwdV[i]);
	   P->mfwdV[i] = L-beg[i];
	   P->BfwdV[i] = malloc (256 * sizeof(Mask));
	   P->SfwdV[i] = malloc (256 * sizeof(Mask));
           extendedLoadVerif (P->mfwdV[i],B,S,A,L-P->mfwdV[i],+1,
	                      P->BfwdV[i],P->SfwdV[i],&P->IfwdV[i],
			      &P->FfwdV[i],&P->AfwdV[i],&P->DifwdV[i]);
	 }
     P->V1 = malloc ((2*(k+1))*sizeof(Mask));
     for (i=0;i<2*(k+1);i++) P->V1[i] = createMask(L);
     P->V2 = createMasks (2,L);

     	/* fifth, the data for fast scanning strongly depends on the
     	   type of search used */

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
     if (K == 0) 
        { P->type = wlen1 ? BWD : FWD; 
	  for (i=beg1;i<end1;i++) SET(active,i);
	}
     else 
	{ P->type = KPLUS1;
	  for (l=0;l<=k;l++)
	     for (i=beg[l];i<end[l];i++) SET(active,i);
	}
     if (detClass(tree,active,L) == EXTENDED)
        { P->scanData = eextendedLoadFast (K?wlen:wlen1,k,P->type,B,S,A,
                                           K?beg:&beg1,K?end:&end1);
          P->scanText = (escanProc)eextendedScan;
          P->scanFree = (freeScanProc)eextendedFreeScan;
        }
     else
        { P->scanData = esimpleLoadFast (K?wlen:end1-beg1,k,P->type,B,
				         K?beg:&beg1);
          P->scanText = (escanProc)esimpleScan;
          P->scanFree = (freeScanProc)esimpleFreeScan;
        }

	/* finally, free temporary structures and return P */

     free (B[0]); free (S[0]); free (A); 
     free (beg); free (end); free (active);
     return P;
   }

void eextendedFree (eextendedData *P)

	/* Frees P */

   { int K,i;
     P->scanFree (P->scanData);
     K = (P->type == KPLUS1) ? P->k : 0;
     for (i=0;i<=K;i++)
         { free (P->BbwdV[i][0]); free (P->SbwdV[i][0]);
           free (P->BbwdV[i]); free (P->SbwdV[i]);
           free (P->IbwdV[i]); free (P->FbwdV[i]); free (P->AbwdV[i]);
	   free (P->DibwdV[i]);
	 }
     free (P->BbwdV); free (P->SbwdV);
     free (P->IbwdV); free (P->FbwdV); free (P->AbwdV);
     free (P->DibwdV);
     for (i=0;i<=K;i++)
         { free (P->BfwdV[i][0]); free (P->SfwdV[i][0]);
           free (P->BfwdV[i]); free (P->SfwdV[i]);
           free (P->IfwdV[i]); free (P->FfwdV[i]); free (P->AfwdV[i]);
	   free (P->DifwdV[i]);
	 }
     free (P->BfwdV); free (P->SfwdV);
     free (P->IfwdV); free (P->FfwdV); free (P->AfwdV);
     free (P->DifwdV);
     for (i=0; i<2*(P->k+1); i++) free (P->V1[i]);
     free (P->V1);
     free (P->V2);
     free (P);
   }

static byte *dirCheck (byte *ptr, byte *top, int inc, 
		       int m, int *K, Mask *B, Mask *S,
		       Mask I, Mask F, Mask A, Mask Di,
		       Mask *D, Mask *T, Mask oD, Mask nD)

		/* directional check, forward or backward.
		   receives maximal error level in *K and returns there
		   the best k found, being the return value the last
		   character considered in that case (NULL if nothing
		   can be found with *K errors) */
		/* assumes k > 0 and that ptr needs to be shifted before
		   reading */
	        /* note that this subpattern CAN have ?s or *s at the
	           beginning (but still not at the end) */

   { mask E,X,Df,nT;		/* masks */ 
     Mask oB,nB,oS,nS;		/* masks */ 
     int i,r,n;
     mask c1,c2,oh1,nh1,oh2,nh2,oh3,nh3;
     byte *ret;
     int k = *K;
     
           /* first remove trivial case m=0 */
      
     if (m == 0)
        { *K = 0;
          while (*K <= k)
	     { if ((inc == -1) && recCheckLeftContext(ptr,top)) return ptr;
	       if ((inc == +1) && recCheckRightContext(ptr+1,top+1)) return ptr;
	       if (ptr == top) return NULL; ptr += inc;
               if (OptIns) (*K)++; else return NULL;
             }
          return NULL;
        }

	  /* multimasks not used for efficiency */

     	   /* fill initial masks, T[i] = 0, D[0] = 1s inherited by ?s
     	      D[i] = D[0] (propagating the 1s if OptDel) */

     ret = NULL;
     n = (m+W-1)/W;
     E = ONE << ((m-1)%W);
     X = ONE;
     for (r=0;r<n;r++) 
	{ T[0][r] = ZEROS; D[0][r] = Di[r]; }
     for (i=1;i<=k;i++)
	{ oh1 = ONE;
	  c1 = ZEROS;
	  for (r=0;r<n;r++) 
	    { T[i][r] = ZEROS; 
	      if (OptDel)
	         { nh1 = (D[i-1][r] >> (W-1)) & ONE;
	           D[i][r] = (D[i-1][r] << 1) | oh1;
		   oh1 = nh1;
	           Df = D[i][r] | F[r];
	           D[i][r] |= A[r] & ((~(Df-c1-I[r])) ^ Df);
	           c1 = (Df < c1 + I[r]) || (c1 + I[r] < c1);
		 }
	      else D[i][r] = D[i-1][r];
	    }
	  if (E & D[i][n-1])  /* k >= m, reasonable here */
	     { if (((inc == -1) && recCheckLeftContext(ptr,top)) ||
	           ((inc == +1) && recCheckRightContext(ptr+1,top+1)))
	          { *K = i; k = i-1; ret = ptr; }
	     } 
	}

     	   /* now process until automaton dies. each time it finds
     	      P with i<=k errors, annote i and pos and reduce k */

     if (ptr == top) return ret;
     ptr += inc; oB = nB = B[*ptr]; oS = nS = S[*ptr];
     while (true)
	{ if (ptr-inc == top) return ret;  /* end of buffer */
	  if (ptr != top) { ptr += inc; nB = B[*ptr]; nS = S[*ptr]; }
	  else ptr += inc;
	     /* oD = D[0]; nD = (((oD << 1)|X) & B[c]) | (oD & S[c]); */
	  oh1 = X;
	  c1 = ZEROS;
	  for (r=0;r<n;r++) 
	      { oD[r] = D[0][r];
		nh1 = (oD[r] >> (W-1)) & ONE;
		nD[r] = (((oD[r] << 1) | oh1) & oB[r]) | (oD[r] & oS[r]);
		oh1 = nh1;
	           /* Df = nD | F; D[0] = nD |= A & ((~(Df-I)) ^ Df); */
	        Df = nD[r] | F[r];
		nD[r] |= A[r] & ((~(Df-c1-I[r])) ^ Df);
		c1 = (Df < c1 + I[r]) || (c1 + I[r] < c1);
		D[0][r] = nD[r];
	      }
	  if (E & nD[n-1])  /* found with 0 errors, return */
	     { if (((inc == -1) && recCheckLeftContext(ptr,top)) ||
	           ((inc == +1) && recCheckRightContext(ptr+1,top+1))) 
		  { *K = 0; return ptr; }
	     }
	  for (i=1;i<=k;i++) 
	      { oh1 = ZEROS;
		oh2 = oh3 = X;
		c1 = c2 = ZEROS;
		for (r=0;r<n;r++) 
		   {     /* nD = nD << 1 (if OptDel, else ZEROS) */
		     if (OptDel)
		        { nh1 = (nD[r] >> (W-1)) & ONE;
	                  nD[r] = (nD[r] << 1) | oh1;
	                  oh1 = nh1;
	                }
		     else nD[r] = ZEROS;
	                 /* nD |= oD (if OptIns) */
		     if (OptIns) nD[r] |= oD[r];
	                 /* nD |= ((oD << 1)|X) (if OptSubs) */
		     nh2 = (oD[r] >> (W-1)) | ONE;
		     if (OptSubs) nD[r] |= (oD[r] << 1) | oh2;
	                 /*nD |= (((D[i] << 1)|X) & B[c])|(D[i] & S[c]) (alws)*/
		     nh3 = (D[i][r] >> (W-1)) & ONE;
	             nD[r] |= (((D[i][r]<<1)|oh3) & oB[r]) | (D[i][r] & oS[r]);
	             oh3 = nh3;
		         /* nD |= T[i] (if OptTransp) */
		     if (OptTransp) 
		         { nD[r] |= T[i][r];
                               /* T[i] = (((oD << 1) | X) & nB) | (oD & nS); */
                           nT = (((oD[r]<<1) | oh2) & nB[r]) | (oD[r] & nS[r]);
			       /* Df = T[i] | F;
				  T[i] |= A & ((~(Df-I)) ^ Df); */
			   Df = nT | F[r];
			   nT |= A[r] & ((~(Df-c2-I[r])) ^ Df);
		           c2 = (Df < c2 + I[r]) || (c2 + I[r] < c2);
			       /* T[i] = ((T[i] << 1) & oB) | (T[i] & oS); */
			   T[i][r] = ((nT << 1) & oB[r]) | (nT & oS[r]);
	                 }
	             oh2 = nh2;
		         /* Df = nD | F; 
		            nD |= (A & ((~(Df-I)) ^ Df)); (always) */
		     Df = nD[r] | F[r];
		     nD[r] |= A[r] & ((~(Df-c1-I[r])) ^ Df);
		     c1 = (Df < c1 + I[r]) || (c1 + I[r] < c1);
		         /* oD = D[i]; D[i] = nD; */
		     oD[r] = D[i][r]; 
		     D[i][r] = nD[r];
		   }
		if (E & nD[n-1])  /* found with k errors, be more strict now */
		   { if (((inc == -1) && recCheckLeftContext(ptr,top)) ||
		         ((inc == +1) && recCheckRightContext(ptr+1,top+1)))
			{ ret = ptr;
		          do i--;
		          while ((i >= 0) && (E & D[i][n-1]));
		          *K = i+1;
		          if (i == -1) return ptr; /* found with 0 errors */
		          k = i;
		        }
		   }
	      }
	  for (r=0;r<n;r++) 
	      { if (r<n-1) { if (D[k][r] | T[k][r]) break; }
		else { if ((D[k][r]  | T[k][r]) & ((E << 1) - ONE)) break; }
	      }
	  if (r == n) return ret;  /* automaton died, return best */
	  X = ZEROS;
	  oB = nB; oS = nS;
        }
   }

static bool checkMatch1 (eextendedData *P,int i, byte *pos, 
		         byte *beg, byte *end)

	/* checkMatch that does not care about a transposition where the
	   pattern is divided */

   { byte *ptrb,*ptrf;		/* ptr to match of selected subpattern */
     int kb,kf;

	/* first part: backward check */

     kb = P->k;
     ptrb = dirCheck (pos,beg,-1,P->mbwdV[i],&kb,P->BbwdV[i],P->SbwdV[i],
		      P->IbwdV[i],P->FbwdV[i],P->AbwdV[i],P->DibwdV[i],
		      P->V1,P->V1+(P->k+1),
		      P->V2,P->V2+(P->mbwdV[i]+P->mfwdV[i]+W-1)/W);
     if (ptrb == NULL) return false;

	/* second part: forward check */

     kf = P->k-kb;
     ptrf = dirCheck (pos-1,end-1,+1,P->mfwdV[i],&kf,P->BfwdV[i],P->SfwdV[i],
		      P->IfwdV[i],P->FfwdV[i],P->AfwdV[i],P->DifwdV[i],
		      P->V1,P->V1+(P->k+1),
		      P->V2,P->V2+(P->mbwdV[i]+P->mfwdV[i]+W-1)/W);
     if (ptrf == NULL) return false;

	/* match found */

     return true;   /* match in [ptrb,ptrf] */
   }

static bool checkMatch (eextendedData *P, int i, byte *pos, 
		        byte **beg, byte **end)

	/* Checks that an initially promising match starting at pos
	   corresponds or not to a complete match. The buffer limits
	   are *beg and *end, and in case of successful match we return
	   there the boundaries of the match found. */

   { byte *rbeg,*rend,*obeg,*oend;
     bool ret;
		
                  /* start by knowing my surrounding record */
      
     if (P->type == FWD)  /* pos-1 is last position read */
        { recGetRecord (pos-1,*beg,*end,&rbeg,&rend,&obeg,&oend);
          if ((rbeg > pos-1) || (rend <= pos-1)) return false; /* overlaps */
        }
     else
        { recGetRecord (pos,*beg,*end,&rbeg,&rend,&obeg,&oend);
          if ((rbeg > pos) || (rend <= pos)) return false; /* overlaps */
        }                                                                       

     ret = checkMatch1(P,i,pos,rbeg,rend);
     if (ret) { *beg = obeg; *end = oend; return true; }

       /* if the pattern is really split and transpositions are permitted
          a delicate problem may appear if the solution needs to transpose
	  the characters at the split point. we solve the problem in a
	  rather brutal but effective way */

     if (OptTransp && (P->mbwdV[i] > 0) && (P->mfwdV[i] > 0) &&
		      (pos < rend) && (pos-1 >= rbeg))
	{ byte c = *pos;
	  *pos = pos[-1];
	  pos[-1] = c;
	  ret = checkMatch1(P,i,pos,rbeg,rend);
          pos[-1] = *pos;
	  *pos = c;
	  if (ret) { *beg = obeg; *end = oend; return true; }
	}

     return false;
   }

static bool kplus1Scan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		        register mask *B, register mask *B0, register mask *S,
			register mask I, register mask F, register mask A,
		        register int m, register int k, mask *ends)

   { register byte *pos,*top;	/* current and final text pointers */
     register mask D,Df;	/* D mask for scanning */
     register int c,j;		/* window index */
     int i;

     pos = *beg; top = *end;

     pos--;
     top -= m;
     while (pos < top)
        { D = B[pos[m]]; if (!D) { pos += m; continue; }  /* skip loop */
	  j = m-1;
          do { Df = D | F;
	       D |= A & ((~(Df-I)) ^ Df);
	       c = pos[j--];
	       D = ((D << 1) & B0[c]) | (D & S[c]); 
	     }
	  while (D && j);
          if (D)
	     {  /* a subpattern has matched, check the rest */
	       for (i=0;i<=k;i++)
		   if ((D & ends[i]) && checkMatch (P,i,pos+1,beg,end))
		      return true;
	 	/* else shift in 1 */
	     }
	  pos += j+1;
        }
     return false;
   }

static bool bwdScan1 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, register mask *S, register mask I, 
		      register mask F, register mask A, int M, register int m)

	/* specialized bwdScan for k=1 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask O,Df;
     register mask D0,D1,T1;
     register int c,j;		/* window index */
     register mask oB,nB,oS,nS;

     pos = *beg; top = *end;
     
     O = ONES << (W-M);
     m -= 2;
     top -= m;
     while (pos < top)
        { j = m;
	     	/* first step is special */
          D0 = B[pos[j--]];
	  c = pos[j];
	  T1 = oB = B[c]; oS = S[c];
	  Df = T1 | F;
	  T1 |= A & ((~(Df-I)) ^ Df);
	  T1 = (T1 << 1) & D0;
	  Df = D0 | F;
	  D0 |= A & ((~(Df-I)) ^ Df);
	  D1 = O;
	 	/* now more steps */
	  while (true)
             { if (j-- > 0) { c = pos[j]; nB = B[c]; nS = S[c]; }
	       D1 = ((D1 << 1) & oB) | (D1 & oS) | D0 | (D0 << 1) | T1;
	       T1 = ((D0 << 1) & nB) | (D0 & nS);
	       Df = T1 | F;
	       T1 |= A & ((~(Df-I)) ^ Df);
	       T1 = (T1 << 1) & oB;
	       D0 = ((D0 << 1) & oB)  | (D0 & oS);
	       Df = D0 | F;
	       D0 |= A & ((~(Df-I)) ^ Df);
	       D1 |= D0 << 1;
	       Df = D1 | F;
	       D1 |= A & ((~(Df-I)) ^ Df);
               if (j < 0) 
	          {  /* selected part of pattern has matched, check rest */
	            if ((D1 & (ONE<<(W-1))) && checkMatch (P,0,pos+1,beg,end)) 
		       return true;
	 	     /* else shift in 1 */
		    break;
	          }
	       if (!D1 && !T1) break;
	       oB = nB; oS = nS;
	     }
	  pos += j+2;
        }
     return false;
  }

static bool fwdScan1 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, register mask *S, register mask I, 
		      register mask F, register mask A, int m)

	/* specialized fwdScan for k=1 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask E,Df;
     register mask D0,D1,T1;
     register mask oB,nB,oS,nS;
     register int c;

     pos = *beg; top = *end;
     if (pos == top) return false;
     
     E = ONE<<(m-1);
     while (pos <= top)
	{ if ((c = *pos++) == OptRecChar) continue;
          D0 = ZEROS; D1 = ONE; T1 = ZEROS;
          oS = S[c]; oB = B[c];
          while (pos <= top)
             { if (pos < top)
		  { if ((c = *pos++) == OptRecChar) break;
		     nS = S[c]; nB = B[c]; 
		  }
	       else pos++;
	       D1 = (((D1 << 1) | ONE) & oB) | (D1 & oS)
	            | D0 | ((D0 << 1) | ONE) | T1;
	       T1 = (((D0 << 1) | ONE) & nB) | (D0 & nS);
	       Df = T1 | F;
	       T1 |= A & ((~(Df-I)) ^ Df);
	       T1 = ((T1 << 1) & oB) | (T1 & oS);
	       D0 = (((D0 << 1) | ONE) & oB)  | (D0 & oS);
	       Df = D0 | F;
	       D0 |= A & ((~(Df-I)) ^ Df);
	       D1 |= D0 << 1;
	       Df = D1 | F;
	       D1 |= A & ((~(Df-I)) ^ Df);
               if (D1 & E)
	          {  /* selected part of pattern has matched, check the rest */
	            if (checkMatch (P,0,pos-1,beg,end)) return true;
	          }
	       oB = nB; oS = nS;
             }
        }
     return false;
   }

static bool bwdScan2 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, register mask *S, register mask I, 
		      register mask F, register mask A, int M, register int m)

	/* specialized bwdScan for k=2 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask O,Df;
     register mask D0,D1,D2,T1,T2;
     register int c,j;		/* window index */
     register mask oB,nB,oS,nS;

     pos = *beg; top = *end;
     
     O = ONES << (W-M);
     m -= 3;
     top -= m;
     while (pos < top)
        { j = m;
     	     /* first step is special */
          D0 = B[pos[j--]];
	  c = pos[j];
	  T1 = oB = B[c]; oS = S[c];
	  Df = T1 | F;
	  T1 |= A & ((~(Df-I)) ^ Df);
	  T1 = (T1 << 1) & D0;
	  Df = D0 | F;
	  D0 |= A & ((~(Df-I)) ^ Df);
	  D1 = D2 = O; T2 = T1;
	     /* now more steps */
	  while (true)
             { if (j-- > 0) { c = pos[j]; nB = B[c]; nS = S[c]; }
	       D2 = ((D2 << 1) & oB) | (D2 & oS) | D1 | (D1 << 1) | T2;
	       T2 = ((D1 << 1) & nB) | (D1 & nS);
	       Df = T2 | F;
	       T2 |= A & ((~(Df-I)) ^ Df);
	       T2 = (T2 << 1) & oB;
	       D1 = ((D1 << 1) & oB) | (D1 & oS) | D0 | (D0 << 1) | T1;
	       T1 = ((D0 << 1) & nB) | (D0 & nS);
	       Df = T1 | F;
	       T1 |= A & ((~(Df-I)) ^ Df);
	       T1 = (T1 << 1) & oB;
	       D0 = ((D0 << 1) & oB) | (D0 & oS);
	       Df = D0 | F;
	       D0 |= A & ((~(Df-I)) ^ Df);
	       D1 |= D0 << 1;
	       Df = D1 | F;
	       D1 |= A & ((~(Df-I)) ^ Df);
	       D2 |= D1 << 1;
	       Df = D2 | F;
	       D2 |= A & ((~(Df-I)) ^ Df);
               if (j < 0) 
	          {  /* selected part of pattern has matched, check rest */
	            if ((D2 & (ONE<<(W-1))) && checkMatch (P,0,pos+1,beg,end)) 
		       return true;
	 	     /* else shift in 1 */
		    break;
	          }
	       if (!D2 && !T2) break;
	       oB = nB; oS = nS;
	     }
	  pos += j+2;
        }
     return false;
   }

static bool fwdScan2 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, register mask *S, register mask I, 
		      register mask F, register mask A, int m)

	/* specialized fwdScan for k=2 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask E,Df;
     register mask D0,D1,D2,T1,T2;
     register mask oB,nB,oS,nS;
     register int c;

     pos = *beg; top = *end;
     if (pos == top) return false;
     
     E = ONE<<(m-1);
     while (pos <= top)
	{ if ((c = *pos++) == OptRecChar) continue;
          D0 = ZEROS; D1 = ONE; D2 = ~(ONES << 2);
          T1 = T2 = ZEROS;
          oS = S[c]; oB = B[c];
          while (pos <= top)
             { if (pos < top)
		  { if ((c = *pos++) == OptRecChar) break;
		     nS = S[c]; nB = B[c]; 
		  }
	       else pos++;
	       D2 = (((D2 << 1) | ONE) & oB) | (D2 & oS) | D1 | (D1 << 1) | T2;
	       T2 = (((D1 << 1) | ONE) & nB) | (D1 & nS);
	       Df = T2 | F;
	       T2 |= A & ((~(Df-I)) ^ Df);
	       T2 = ((T2 << 1) & oB) | (T2 & oS);
	       D1 = (((D1 << 1) | ONE) & oB) | (D1 & oS)
	            | D0 | ((D0 << 1) | ONE) | T1;
	       T1 = (((D0 << 1) | ONE) & nB) | (D0 & nS);
	       Df = T1 | F;
	       T1 |= A & ((~(Df-I)) ^ Df);
	       T1 = ((T1 << 1) & oB) | (T1 & oS);
	       D0 = (((D0 << 1) | ONE) & oB) | (D0 & oS);
	       Df = D0 | F;
	       D0 |= A & ((~(Df-I)) ^ Df);
	       D1 |= D0 << 1;
	       Df = D1 | F;
	       D1 |= A & ((~(Df-I)) ^ Df);
	       D2 |= D1 << 1;
	       Df = D2 | F;
	       D2 |= A & ((~(Df-I)) ^ Df);
               if (D2 & E)
	          {  /* selected part of pattern has matched, check the rest */
	            if (checkMatch (P,0,pos-1,beg,end)) return true;
	          }
	       oB = nB; oS = nS;
             }
        }
     return false;
   }

static bool bwdScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		     register mask *B, register mask *S, register mask I, 
		     register mask F, register mask A, register mask *D, 
		     register mask *T, int M, register int m, register int k)

	   /* backward scanning with suffix automaton */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask oD,nD,nT;	/* old, new, current D mask */
     register mask O,Df;		/* masks */
     register int c,i,j;		/* window index */
     register mask oB,nB,oS,nS;

     pos = *beg; top = *end;
     
     O = ONES << (W-M);
     m -= k+1;
     top -= m;
     while (pos < top)
        { j = m;
	     	/* first step is special: there is a dirty play with
	     	   names to save some reassignments */
          nD = B[pos[j--]];
	  c = pos[j];
	  nT = oB = B[c]; oS = S[c];
	  Df = nT | F;
	  nT |= A & ((~(Df-I)) ^ Df);
	  nT = (nT << 1) & nD;
	  Df = nD | F;
	  nD |= A & ((~(Df-I)) ^ Df);
	  D[0] = nD;
          for (i=1;i<=k;i++) { D[i] = O; T[i] = nT; }
       		/* now more steps */
          while (true)
            { if (j-- > 0) { c = pos[j]; nB = B[c]; nS = S[c]; }
	      oD = D[0]; 
	      nD = ((oD << 1) & oB) | (oD & oS);
	      Df = nD | F;
	      nD |= A & ((~(Df-I)) ^ Df);
	      D[0] = nD;
	      for (i=1;i<=k;i++)
	          { nD = ((D[i] << 1) & oB) | (D[i] & oS)
	                 | oD | ((oD | nD) << 1) | T[i];
		    Df = nD | F;
		    nD |= A & ((~(Df-I)) ^ Df);
		    nT = ((oD << 1) & nB) | (oD & nS);
		    Df = nT | F;
		    nT |= A & ((~(Df-I)) ^ Df);
		    T[i] = (nT << 1) & oB;
		    oD = D[i]; D[i] = nD;
		  }
              if (j < 0) 
	         {  /* selected part of pattern has matched, check rest */
	           if ((nD & (ONE<<(W-1))) && checkMatch (P,0,pos+1,beg,end)) 
		      return true;
	 	    /* else shift in 1 */
		   break;
	         }
	      if (!nD && !T[k]) break;
	      oB = nB; oS = nS;
	    }
	  pos += j+2;
        }
     return false;
   }

static bool fwdScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		     register mask *B, register mask *S, register mask I, 
		     register mask F, register mask A, register mask *D, 
		     register mask *T, int m, register int k)

	   /* forward scanning */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask oD,nD,nT;	/* old, new D,T masks */
     register mask Df;		/* masks */
     register int i,c;		/* window index */
     register mask E,oB,nB,oS,nS;

     pos = *beg; top = *end;
     if (pos == top) return false;
     
     E = ONE<<(m-1);
     while (pos <= top)
	{ if ((c = *pos++) == OptRecChar) continue;
          for (i=0;i<=k;i++) { D[i] = ~(ONES << i); T[i] = ZEROS; }
          oS = S[c]; oB = B[c];
          while (pos <= top)
             { if (pos < top) 
		  { if ((c = *pos++) == OptRecChar) break;
		     nS = S[c]; nB = B[c]; 
		  }
	       else pos++;
	       oD = D[0];
	       nD = ((oD << 1) | ONE) & oB;
	       Df = nD | F;
	       nD |= A & ((~(Df-I)) ^ Df);
	       D[0] = nD;
	       for (i=1;i<=k;i++)
	           { nD = (((D[i] << 1) | ONE) & oB) | (D[i] & oS)
	                  | oD | (((oD | nD) << 1) | ONE) | T[i];
		     Df = nD | F;
		     nD |= A & ((~(Df-I)) ^ Df);
		     nT = (((oD << 1) | ONE) & nB) | (oD & nS);
		     Df = nT | F;
		     nT |= A & ((~(Df-I)) ^ Df);
		     T[i] = ((nT << 1) & oB) | (nT & oS);
		     oD = D[i]; D[i] = nD;
	           }
               if (nD & E)
	          {  /* selected part of pattern has matched, check the rest */
	            if (checkMatch (P,0,pos-1,beg,end)) return true;
	          }
	       oB = nB; oS = nS;
             }
        }
     return false;
   }

bool eextendedScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		    eextendedScanData *scan)

	/* Scans from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P. It requires
	     - the type of search "type" to perform
	     - a procedure 
	        bool checkMatch (P,byte *pos, byte **beg, byte **end)
	       which whether there is a match around pos and returns in
               *beg and *end the limits. "Around" means that pos is the initial
	       position of the search pattern given to esimpleScan.
	     - tables of masks B,S,I,F,A for the pattern
	     - for KPLUS1, a table of masks B0, =B but killing initial positions
	     - for !KPLUS1, space for k+1 masks in D and T
	     - for KPLUS a table of ends masks
	     - length M of the search pattern(s) and number k of errors
	     - window length m
	 */

   { switch (scan->type)
        { case KPLUS1: 
	    return kplus1Scan (beg,end,checkMatch,P,scan->edata->B,
				scan->B0,scan->edata->S,scan->edata->I,
				scan->edata->F,scan->edata->A,
				scan->edata->wlen,scan->k,scan->ends);
	  case FWD: 
	    if (scan->k==1) return fwdScan1(beg,end,checkMatch,P,scan->edata->B,
					 scan->edata->S,scan->edata->I,
					 scan->edata->F,scan->edata->A,scan->m);
	    if (scan->k==2) return fwdScan2(beg,end,checkMatch,P,scan->edata->B,
					 scan->edata->S,scan->edata->I,
					 scan->edata->F,scan->edata->A,scan->m);
	    return fwdScan (beg,end,checkMatch,P,scan->edata->B,scan->edata->S,
			     scan->edata->I,scan->edata->F,scan->edata->A,
			     scan->V3,scan->V3+(scan->k+1),scan->m,scan->k);
	  case BWD: 
	    if (scan->k==1) return bwdScan1(beg,end,checkMatch,P,scan->edata->B,
					 scan->edata->S,scan->edata->I,
					 scan->edata->F,scan->edata->A,scan->m,
					 scan->edata->wlen);
	    if (scan->k==2) return bwdScan2(beg,end,checkMatch,P,scan->edata->B,
					 scan->edata->S,scan->edata->I,
					 scan->edata->F,scan->edata->A,scan->m,
					 scan->edata->wlen);
	    return bwdScan (beg,end,checkMatch,P,scan->edata->B,scan->edata->S,
			     scan->edata->I,scan->edata->F,scan->edata->A,
			     scan->V3,scan->V3+(scan->k+1),scan->m,
			     scan->edata->wlen,scan->k);
        }
   }

bool eextendedSearch (byte **beg, byte **end, eextendedData *P)

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

   { return P->scanText (beg,end,(echeckProc)checkMatch,P,P->scanData);
   }

