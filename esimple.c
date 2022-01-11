
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

	/* Search a simple pattern in a text allowing errors */

#include "esimple.h"
#include "parser.h"
#include "options.h"

#define BK 0.97  /* k+1 has to be at most this to prefer it over fwd/bwd */

static double findBestMulti (Mask *B, int m, int K, int *wlen, int *beg)

	/* Finds a best set of K+1 factors to search, puts their length
	   in *wlen and their initial positions in beg[i]. If no good
	   set exists it returns wlen=0
	*/

   { int tr = OptTransp ? 1 : 0;
     int L = min(m-tr*K,W)/(K+1);
     double *prob;  /* prob[i] = prob pattern letter i */
     double *pprob; /* pprob[i,l] = matching prob of a factor
		                    of length l starting at position i */
     double *pcost; /* pcost[i,l] = average cost per window of a factor
		                    of length l starting at i */
     double *mprob; /* for fixed i, mprob[l,k] = prob of any factor of
		       length k in P[i..i+l-1] */
     double *mbest; /* mbest[i,k] = cost using best selection of k strings 
		       starting at i. the length of the strings is fixed */
     int *ibest;    /* where to start the first strings to obtain the
		       mbest cost */
     int c,i,k,l,wl;
     double this,best,cost;

		/* first, load probabilities, O(m) time */

     prob = malloc (m * sizeof (double));
     for (i=0; i<m; i++) 
       { prob[i] = 0.0;
	 for (c=0; c<256; c++)
	     if (ISSET(B[c],i)) prob[i] += letterProb[c];
       }

     		/* second, compute pprob[i,l], O(m^2/k) time */

     pprob = malloc ((m+1) * (L+1) * sizeof (double));
     pprob[(L+1)*m+0] = 1.0;
     for (l=1;l<=L;l++) pprob[(L+1)*m+l] = 0.0;
     for (i=m-1;i>=0;i--)
         { pprob[(L+1)*i+0] = 1.0;
	   for (l=1;l<=L;l++)
               pprob[(L+1)*i+l] = prob[i] * pprob[(L+1)*(i+1)+(l-1)];
	 }
     free (prob);
	
     		/* third, compute pcost[i,l], O(m^3/k^2) time */
     		   
     mprob = malloc ((L+1) * L * sizeof (double));
     pcost = malloc (m * L * sizeof (double));
     for (i=0; i<m; i++)
	{ for (k=1;k<=L;k++) mprob[(0)*L+(k-1)] = 0.0;
	  for (l=1; l<=L; l++)
	     { cost = 1.0;
	       for (k=1;k<=l;k++)  /* length of scanned window */
		  { mprob[l*L+(k-1)] = 1-
		       (1-mprob[(l-1)*L+(k-1)])*(1-pprob[(L+1)*(i+l-k)+k]);
		    cost += mprob[l*L+(k-1)];
	          }
	       pcost[i*L+(l-1)] = cost;
	     }
	}
     free (pprob);
     free (mprob);

		/* fourth, find best factors. this takes O(m^2) time in
		   the worst case but the average should be closer to O(km) */
    
     mbest = malloc ((m+1) * (K+2) * sizeof (double));
     ibest = malloc ((m+1) * (K+2) * sizeof (int));
     *wlen = 0; best = BK;
     for (wl = L; wl > 1; wl--)
	{ if (best < 1/(double)wl) break;  /* cannot win */
	  for (i=0;i<=m;i++) mbest[i*(K+2)+(0)] = 0.0;
	  for (k=1;k<=K+1;k++) mbest[m*(K+2)+k] = 1.0;
	  for (k=1; k<=K+1; k++)
	     for (i=m-k*wl-(k-1)*tr; i>=0; i--)
                 { cost = (pcost[i*L+(wl-1)] < wl+1) ?
			  pcost[i*L+(wl-1)]/(wl-pcost[i*L+(wl-1)]+1) : 1.0;
		   if (cost > 1.0) cost = 1.0;
		   this = 1-(1-cost)*(1-mbest[(i+wl+tr)*(K+2)+(k-1)]);
		   ibest[i*(K+2)+k] = i;
	           if ((i+1<=m-k*wl-(k-1)*tr) && (mbest[(i+1)*(K+2)+k] < this))
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
		    i = beg[K+1-k] + wl + tr;
		  }
	       best = mbest[0*(K+2)+(K+1)];
	     }
	}
     free (pcost);
     free (mbest);
     free (ibest);

     if (best < BK) return best;
     *wlen = 0;
     return 1.0;
   }

void esimpleFreeScan (esimpleScanData *scan)

   { simpleFreeScan(scan->sdata);
     free (scan->V3);
     free (scan);
   }

esimpleScanData *esimpleLoadFast (int wlen, int k, int type, Mask *B, int *beg)

   { int i,c,l;
     esimpleScanData *scan = malloc (sizeof(esimpleScanData));
     scan->k = k;
     scan->type = type;
     if (type != KPLUS1)
        { scan->sdata = simpleLoadFast (type == BWD,B,*beg,*beg+wlen);
          scan->V3 = malloc ((2*(k+1))*sizeof(mask));
        }
     else   /* partition in k+1 */
        { scan->sdata = malloc (sizeof(simpleScanData));
          scan->sdata->m = wlen;
          scan->V3 = NULL;
          for (c=0; c<256; c++) scan->sdata->B[c] = ZEROS;
          for (i=0;i<=k;i++)
              { for (l=0; l<wlen; l++)
	           for (c=0; c<256; c++) 
	              if (ISSET(B[c],beg[i]+wlen-1-l))
		         { scan->sdata->B[c] |= ONE << (wlen*i+l); 
		           if (l > 0) scan->B0[c] |= ONE << (wlen*i+l); 
		         }
              }
        }
     return scan;
   }

esimpleData *esimplePreproc (byte *pat, Tree *tree, Tree **pos, int k)

	/* Preprocesses pat and creates an esimpleData structure
	   for searching */

   { esimpleData *P = malloc (sizeof(esimpleData));
     Mask B[256];
     int r,i,c,l,L,wlen;
     bool bwd;
     int *beg,beg1,end1;
     double best1,best2;
     int K;

        /* first, compute the length of the pattern in bits */ 

     L = 0;
     simpleLength (tree,&L);

        /* second, create bit masks */

     simpleLoadMasks (pat,strlen(pat),L,B,pos);

	/* third, find the best way to search the pattern */

     best1 = (k+1) * simpleFindBest (B,L,k,&bwd,&beg1,&end1); 
     beg = malloc ((k+1)*sizeof(int));
     best2 = findBestMulti (B,L,k,&wlen,beg); 
     if ((best1 <= best2) || !wlen)  /* we prefer suffix */
	{ K = 0; beg[0] = bwd ? beg1 : end1; }
     else K = k;
	
	/* fourth, create the verification data, common proc */

     P->k = k;
     P->mbwdV = malloc ((K+1)*sizeof(int));
     P->BbwdV = malloc ((K+1)*sizeof(Mask*));
     P->mfwdV = malloc ((K+1)*sizeof(int));
     P->BfwdV = malloc ((K+1)*sizeof(Mask*));
     for (i=0;i<=K;i++) 
	 { P->mbwdV[i] = beg[i];
	   P->BbwdV[i] = malloc (256 * sizeof(Mask));
	   simpleLoadVerif (P->mbwdV[i],B,P->mbwdV[i]-1,-1,P->BbwdV[i]);
	   P->mfwdV[i] = L-beg[i];
	   P->BfwdV[i] = malloc (256 * sizeof(Mask));
	   simpleLoadVerif (P->mfwdV[i],B,L-P->mfwdV[i],+1,P->BfwdV[i]);
	 }
     P->V1 = malloc ((2*(k+1))*sizeof(Mask));
     for (i=0;i<2*(k+1);i++) P->V1[i] = createMask (L);
     P->V2 = createMasks (2,L);

     	/* fifth, the data for fast scanning strongly depends on the
     	   type of search used */

     if (K == 0) P->type = bwd ? BWD : FWD;
     else P->type = KPLUS1;  /* partition in k+1 */
     P->scanData = esimpleLoadFast (K?wlen:end1-beg1,k,P->type,B,K?beg:&beg1);
     P->scanText = (escanProc)esimpleScan;
     P->scanFree = (freeScanProc)esimpleFreeScan;

	/* finally, free temporary structures and return P */

     free (B[0]);
     free (beg);
     return P;
   }

void esimpleFree (esimpleData *P)

	/* Frees P */

   { int K,i;
     P->scanFree (P->scanData);
     K = (P->type == KPLUS1) ? P->k : 0;
     for (i=0;i<=K;i++)
         { free (P->BbwdV[i][0]);
           free (P->BbwdV[i]);
	 }
     free (P->BbwdV);
     for (i=0;i<=K;i++)
         { free (P->BfwdV[i][0]);
           free (P->BfwdV[i]);
	 }
     free (P->BfwdV);
     for (i=0; i<2*(P->k+1); i++) free (P->V1[i]);
     free (P->V1);
     free (P->V2);
     free (P);
   }

static byte *dirCheck (byte *ptr, byte *top, int inc, 
		       int m, int *K, Mask *B, 
		       Mask *D, Mask *T, Mask oD, Mask nD)

		/* directional check, forward or backward.
		   receives maximal error level in *k and returns there
		   the best k found, being the return value the last
		   character considered in that case (NULL if nothing
		   can be found with *k errors) */
		/* assumes m,k > 0 and that ptr needs to be shifted before
		   reading */

   { mask E,X;			/* masks */ 
     Mask nB;
     int i,r,n;
     mask oh1,nh1,oh2,nh2,oh3,nh3,oh4,nh4,oh5,nh5;
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

     	   /* fill initial masks, T[i] = 0, D[i] = ~(ONES<<i) (if OptDel) */

     ret = NULL;
     n = (m+W-1)/W;
     E = ONE << ((m-1)%W);
     X = ONE;
     for (i=0;i<=k;i++)
	 { for (r=0;r<n;r++) 
	      { T[i][r] = ZEROS; 
	        if (OptDel)
	           D[i][r] = (r+1)*W <= i ? ONES :
		             (r*W >= i ? ZEROS : ~(ONES << (i%W)));
	        else D[i][r] = ZEROS;
	      }
	   if (E & D[i][n-1])  /* k >= m, reasonable here */
	      { if (((inc == -1) && recCheckLeftContext(ptr,top)) ||
	            ((inc == +1) && recCheckRightContext(ptr+1,top+1)))
	           { *K = i; k = i-1; ret = ptr; }
	      }
	 }

     	   /* now process until automaton dies. each time it finds
     	      P with i<=k errors, record i and pos and reduce k */

     while (true)
	{ if (ptr == top) return ret;  /* end of buffer */
	  ptr += inc;
	  nB = B[*ptr]; /* take new character */
	     /* oD = D[0]; nD = D[0] = ((oD << 1)|X) & B[c]; */
	  oh1 = X;
	  for (r=0;r<n;r++) 
	      { oD[r] = D[0][r];
		nh1 = (oD[r] >> (W-1)) & ONE;
		nD[r] = D[0][r] = ((oD[r] << 1) | oh1) & nB[r];
		oh1 = nh1;
	      }
	  if (E & nD[n-1])  /* found with 0 errors, check context & return */
	     { if (((inc == -1) && recCheckLeftContext(ptr,top)) ||
	           ((inc == +1) && recCheckRightContext(ptr+1,top+1)))
	          { *K = 0; return ptr; }
	     }
	  for (i=1;i<=k;i++) 
	      { oh1 = oh4 = ZEROS;
		oh2 = oh3 = X;
		oh5 = X << 1;
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
		     if (OptSubs)
		        { nh2 = (oD[r] >> (W-1)) & ONE;
	                  nD[r] |= (oD[r] << 1) | oh2;
	                  oh2 = nh2;
	                }
	                 /* nD |= ((D[i] << 1)|X) & B[c] (always) */
		     nh3 = (D[i][r] >> (W-1)) & ONE;
	             nD[r] |= ((D[i][r] << 1) | oh3) & nB[r];
	             oh3 = nh3;
		         /* nD |= T[i] & (B[c]<<1); (if OptTransp) */
		     if (OptTransp)
		         { nh4 = (nB[r] >> (W-1)) & ONE;
	                   nD[r] |= ((nB[r] << 1) | oh4) & T[i][r];
	                   oh4 = nh4;
		              /* T[i] = ((oD << 2)|(ONE<<1)) & B[c]; */
		           nh5 = (oD[r] >> (W-2)) & ~(ONES<<2);
	                   T[i][r] = ((oD[r] << 2) | oh5) & nB[r];
	                   oh5 = nh5;
	                 }
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
        }
   }

static bool checkMatch1(esimpleData *P,int i, byte *pos, byte *beg, byte *end)

	/* checkMatch that does not care about a transposition where the
	   pattern is divided */

   { byte *ptrb,*ptrf;		/* ptr to match of selected subpattern */
     int kb,kf;
       
	/* first part: backward check */

     kb = P->k;
     ptrb = dirCheck (pos,beg,-1,P->mbwdV[i],&kb,P->BbwdV[i],
		      P->V1,P->V1+(P->k+1),
		      P->V2,P->V2+(P->mbwdV[i]+P->mfwdV[i]+W-1)/W);
     if (ptrb == NULL) return false;

	/* second part: forward check */

     kf = P->k-kb;
     ptrf = dirCheck (pos-1,end-1,+1,P->mfwdV[i],&kf,P->BfwdV[i],
		      P->V1,P->V1+(P->k+1),
		      P->V2,P->V2+(P->mbwdV[i]+P->mfwdV[i]+W-1)/W);
     if (ptrf == NULL) return false;

	/* match found */

     return true;  /* match in [ptrb,ptrf] */
   }

static bool checkMatch(esimpleData *P, int i, byte *pos, byte **beg, byte **end)

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
		        register mask *B, register mask *B0,
		        register int m, register int k)

   { register byte *pos,*top;	/* current and final text pointers */
     register mask D;		/* D mask for scanning */
     register int j;		/* window index */
     int i;

     pos = *beg; top = *end;

     pos--;
     top -= m;
     while (pos < top)
        { D = B[pos[m]]; if (!D) { pos += m; continue; }  /* skip loop */
	  j = m-1;
          do D = (D << 1) & B0[pos[j--]]; while (D && j);
          if (D)
	     {  /* a subpattern has matched, check the rest */
	       for (i=0;i<=k;i++)
		   if ((D & (1<<(m*(i+1)-1))) && checkMatch (P,i,pos+1,beg,end))
		      return true;
	 	/* else shift in 1 */
	     }
	  pos += j+1;
        }
     return false;
   }

static bool bwdScan1 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, register int m)

	/* specialized bwdScan for k=1 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask O;
     register mask D0,D1,T1;
     register int j;		/* window index */
     register mask B1;

     pos = *beg; top = *end;
     
     O = ONES << (W-m);
     m -= 2;
     top -= m;
     while (pos < top)
        { j = m;
	     	/* first step is special */
	  D0 = T1 = B[pos[j]];
	  D1 = O;
	 	/* now more steps */
	  while (true)
             { B1 = B[pos[--j]];
	       D1 = ((D1 << 1) & B1) | D0 | (D0 << 1) | (T1 & (B1 << 1));
	       T1 = (D0 << 2) & B1;
	       D0 = (D0 << 1) & B1;
	       D1 |= D0 << 1;
               if (j == 0) 
	          {  /* selected part of pattern has matched, check rest */
	            if ((D1 & (ONE<<(W-1))) && checkMatch (P,0,pos,beg,end)) 
		       return true;
	 	     /* else shift in 1 */
		    break;
	          }
	       if (!D1 && !T1) break;
	     }
	  pos += j+1;
        }
     return false;
  }

static bool fwdScan1 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, int m)

	/* specialized fwdScan for k=1 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask E;
     register mask D0,D1,T1;
     register mask B1;

     pos = *beg; top = *end;
     
     E = ONE<<(m-1);
     D0 = ONES; D1 = ONES << 1; T1 = ONES;
     while (pos < top)
        { B1 = B[*pos++];
	  D1 = ((D1<<1) | B1) & D0 & (D0 << 1) & (T1 | ((B1 << 1) | ONE));
	  T1 = (D0 << 2) | B1;
	  D0 = (D0 << 1) | B1;
	  D1 &= D0 << 1;
          if (!(D1 & E))
	     {  /* selected part of pattern has matched, check the rest */
	       if (checkMatch (P,0,pos,beg,end)) return true;
	     }
        }
     return false;
   }

static bool bwdScan2 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, register int m)

	/* specialized bwdScan for k=2 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask O;
     register mask D0,D1,D2,T1,T2;
     register int j;		/* window index */
     register mask B1,B2;

     pos = *beg; top = *end;
     
     O = ONES << (W-m);
     m -= 3;
     top -= m;
     while (pos < top)
        { j = m;
     	     /* first step is special */
	  D0 = T1 = T2 = B[pos[j]];
	  D1 = D2 = O;
	     /* now more steps */
	  while (true)
             { B1 = B[pos[--j]];
	       B2 = B1 << 1;
	       D2 = ((D2 << 1) & B1) | D1 | (D1 << 1) | (T2 & B2);
	       T2 = (D1 << 2) & B1;
	       D1 = ((D1 << 1) & B1) | D0 | (D0 << 1) | (T1 & B2);
	       T1 = (D0 << 2) & B1;
	       D0 = (D0 << 1) & B1;
	       D1 |= D0 << 1;
	       D2 |= D1 << 1;
               if (j == 0) 
	          {  /* selected part of pattern has matched, check rest */
	            if ((D2 & (ONE<<(W-1))) && checkMatch (P,0,pos,beg,end)) 
		       return true;
	 	     /* else shift in 1 */
		    break;
	          }
	       if (!D2 && !T2) break;
	     }
	  pos += j+1;
        }
     return false;
   }

static bool fwdScan2 (byte **beg, byte **end, echeckProc checkMatch, void *P,
		      register mask *B, int m)

	/* specialized fwdScan for k=2 */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask E;
     register mask D0,D1,D2,T1,T2;
     register mask B1,B2;

     pos = *beg; top = *end;
     
     E = ONE<<(m-1);
     D0 = ONES; D1 = ONES << 1; D2 = ONES << 2;
     T1 = T2 = ONES;
     while (pos < top)
        { B1 = B[*pos++];
          B2 = (B1 << 1) | ONE;
	  D2 = ((D2 << 1) | B1) & D1 & (D1 << 1) & (T2 | B2);
	  T2 = (D1 << 2) | B1;
	  D1 = ((D1 << 1) | B1) & D0 & (D0 << 1) & (T1 | B2);
	  T1 = (D0 << 2) | B1;
	  D0 = (D0 << 1) | B1;
	  D1 &= D0 << 1;
	  D2 &= D1 << 1;
          if (!(D2 & E))
	     {  /* selected part of pattern has matched, check the rest */
	       if (checkMatch (P,0,pos,beg,end)) return true;
	     }
        }
     return false;
   }

static bool bwdScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		     register mask *B, register mask *D, register mask *T,
		     register int m, register int k)

	   /* backward scanning with suffix automaton */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask oD,nD;	/* old, new D mask */
     register mask O;		/* masks */
     register int i,j;		/* window index */
     register mask B1,B2;

     pos = *beg; top = *end;
     
     O = ONES << (W-m);
     m -= k+1;
     top -= m;
     while (pos < top)
        { j = m;
	     	/* first step is special */
          D[0] = B[pos[j]];
          for (i=1;i<=k;i++) { D[i] = O; T[i] = D[0]; }
       		/* now more steps */
          while (true)
            { B1 = B[pos[--j]];
	      B2 = B1 << 1;
	      oD = D[0]; D[0] = nD = (oD << 1) & B1;
	      for (i=1;i<=k;i++)
	          { nD = ((D[i]<<1) & B1) | oD | ((oD | nD)<<1) | (T[i] & B2);
		    T[i] = (oD << 2) & B1;
		    oD = D[i]; D[i] = nD;
		  }
              if (j == 0) 
	         {  /* selected part of pattern has matched, check rest */
	           if ((nD & (ONE<<(W-1))) && checkMatch (P,0,pos,beg,end)) 
		      return true;
	 	    /* else shift in 1 */
		   break;
	         }
	      if (!nD && !T[k]) break;
	    }
	  pos += j+1;
        }
     return false;
   }

static bool fwdScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		     register mask *B, register mask *D, register mask *T,
		     int m, register int k)

	   /* forward scanning */

   { register byte *pos,*top;	/* current and final text pointers */
     register mask oD,nD;	/* old, new D mask */
     register mask E;		/* masks */
     register int i,j;	/* window index */
     register mask B1,B2;

     pos = *beg; top = *end;
     
     E = ONE<<(m-1);
     for (i=0;i<=k;i++) { D[i] = ONES << i; T[i] = ONES; }
     while (pos < top)
        { B1 = B[*pos++];
          B2 = (B1 << 1) | ONE;
	  oD = D[0];
	  nD = D[0] = (oD << 1) | B1;
	  for (i=1;i<=k;i++)
	     { nD = ((D[i] << 1) | B1) & oD & ((oD & nD) << 1) & (T[i] | B2);
	       T[i] = (oD << 2) | B1;
	       oD = D[i]; D[i] = nD;
	     }
          if (!(nD & E))
	     {  /* selected part of pattern has matched, check the rest */
	       if (checkMatch (P,0,pos,beg,end)) return true;
	     }
        }
     return false;
   }

bool esimpleScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		  esimpleScanData *scan)

	/* Scans from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P. It requires
	     - the type of search "type" to perform
	     - a procedure 
	        bool checkMatch (P,byte *pos, byte **beg, byte **end)
	       which whether there is a match around pos and returns in
               *beg and *end the limits. "Around" means that pos is the initial
	       position of the search pattern given to esimpleScan.
	     - a table of masks B for the pattern
	     - for KPLUS1, a table of masks B0, =B but killing initial positions
	     - for !KPLUS1, space for k+1 masks in D and T
	     - length m of the search pattern(s) and number k of errors
	 */

   { switch (scan->type)
        { case KPLUS1: 
	    return kplus1Scan (beg,end,checkMatch,P,scan->sdata->B,
			       scan->B0,scan->sdata->m,scan->k);
	  case FWD: 
	    if (scan->k==1) return fwdScan1(beg,end,checkMatch,P,scan->sdata->B,
					    scan->sdata->m);
	    if (scan->k==2) return fwdScan2(beg,end,checkMatch,P,scan->sdata->B,
					    scan->sdata->m);
	    return fwdScan (beg,end,checkMatch,P,scan->sdata->B,scan->V3,
			    scan->V3+(scan->k+1),scan->sdata->m,scan->k);
	  case BWD: 
	    if (scan->k==1) return bwdScan1(beg,end,checkMatch,P,scan->sdata->B,
					    scan->sdata->m);
	    if (scan->k==2) return bwdScan2(beg,end,checkMatch,P,scan->sdata->B,
					    scan->sdata->m);
	    return bwdScan (beg,end,checkMatch,P,scan->sdata->B,scan->V3,
			    scan->V3+(scan->k+1),scan->sdata->m,scan->k);
        }
   }

bool esimpleSearch (byte **beg, byte **end, esimpleData *P)

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

   { return P->scanText (beg,end,(echeckProc)checkMatch,P,P->scanData);
   }

