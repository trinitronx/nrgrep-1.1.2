
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

	/* Search a regular expression in a text allowing errors */

#include "eregular.h"
#include "esimple.h"
#include "eextended.h"
#include "parser.h"

#define BK 0.78  /* k+1 has to be at most this to prefer it over fwd/bwd */

static double findBestMulti (Tree *e, Mask *trans, Mask *rtrans, int m, 
                        Mask *B, Mask initial, Mask final, int *wlen,
                        Mask *dactive, Mask *dinitial, Mask *dfinal, int K)

        /* Finds K+1 best subexpressions to search and puts their states in
	   dactive,dinitial,dfinal and in wlen the window length (0 if it 
	   recommends forward scanning). 
           Returns the avg number of chars inspected. All those masks arrays
	   are allocated, but their contents are not */

   { double *prob;    /* prob[i] = prob of arriving to state i from prev */
     double *pprob;   /* pprob[i,l] = prob matching a path of length l from
                                     state i */
     Mask *initials;  /* initial states at distance i from start state */
     double *mprob;   /* mprob[i,l,l'] = prob of any factor of length l'
                         inside a path of length l departing from state i */
     double *pcost;    /* pcost[i,l] = avg cost per window in same area */
     double *mbest; /* mbest[i,k] = cost using best selection of k
                       strings starting at initials[i]. the length of the 
		       strings is fixed */
     int *ibest;    /* where to start the first strings to obtain the
                       mbest cost */
     Mask reach,tmp;
     int k,i,c,j,wl,l,lp,L,L2;
     int tr = OptTransp ? 1 : 0;
     double this,best,cost;

                /* compute max wlen wl <= m, in O(wl*m^2/w) = 
		   O(m^3/w) time */

     reach = createMask(m);
     tmp = createMask(m);
     wl = 1; /* there is an extra initial state */
     SET(ZERO(reach,m),0);
     while (ISZERO(AND(COPY(tmp,reach,m),final,m),m))
        { wl++;
          ZERO(tmp,m);
          for (j=0;j<m;j++)
             if (ISSET(reach,j)) OR (tmp,trans[j],m);
                OR(reach,tmp,m);
        }
     L = 1+min(wl-1-tr*K,W)/(K+1);  /* 1+max length of pieces */
     if (L <= 1) { *wlen = 0; return 1.0; }
     L2 = L*L;

		/* compute initial states, O(wl*m^2/w) time */

     initials = malloc (wl * sizeof(Mask));
     initials[0] = createMasks (wl+1,m);
     for (i=1;i<wl;i++) initials[i] = initials[0] + i*maskSize(m);
     COPY(initials[0],initial,m); COPY(reach,initial,m);
     for (i=1;i<wl;i++)
        { ZERO(tmp,m);
          for (j=0;j<m;j++)
              if (ISSET(reach,j)) OR(tmp,trans[j],m);
	  AND(NOT(COPY(initials[i],reach,m),m),tmp,m);
	  OR(reach,tmp,m);
        }
     free (tmp); free (reach);
     
                /* load letter probabilities in O(m*s) time */

     prob = malloc (m * sizeof (double));
     for (i=0;i<m;i++)
         { prob[i] = 0.0;
           for (c=0;c<256;c++)
               if (ISSET(B[c],i)) prob[i] += letterProb[c];
         }

                /* compute path probabilities in O(wl*m^2*s/k) time */
                /* pprob[i,l] includes already the prob. of reaching i
                   from its predecessor */

     pprob = malloc (m * L * sizeof (double));
     for (i=0;i<m;i++)
        { pprob[i*L+0] = 1.0;
          pprob[i*L+1] = prob[i];
        }
     for (l=2;l<L;l++)
        for (i=0;i<m;i++)
           { pprob[i*L+l] = 0.0;
             if (!ISZERO(trans[i],m))
                for (j=0;j<m;j++)
                   if (ISSET(trans[i],j))
                      pprob[i*L+l] += pprob[j*L+(l-1)];
             pprob[i*L+l] *= prob[i];
	     if (pprob[i*L+l] > 1.0) pprob[i*L+l] = 1.0;
           }
     free (prob);

                /* compute mprob in O((m*wl/k)^2) time */

     mprob = malloc (m * L2 * sizeof (double));
     for (i=0;i<m;i++)
        for (l=0;l<L;l++)
           mprob[i*L2+l*L+0] = 1.0;
     for (l=1;l<L;l++)
        for (lp=1;lp<=l;lp++)
           for (i=0;i<m;i++)
              { mprob[i*L2+l*L+lp] = pprob[i*L+lp];
                if (lp < l)
                   if (!ISZERO(trans[i],m))
                      for (j=0;j<m;j++)
                         if (ISSET(trans[i],j))
			    mprob[i*L2+l*L+lp] = 1-(1-mprob[i*L2+l*L+lp])*
                                		   (1-mprob[j*L2+(l-1)*L+lp]);
              }
     free (pprob);

                /* compute pcost in O(m*(wl/k)^2) time */

     pcost = malloc (wl * L * sizeof (double));
     for (i=0;i<wl;i++)
        for (l=0;l<L;l++)
           { double new = 1.0;
	     for (j=0;j<m;j++)
		if (ISSET(initials[i],j))
                   for (lp=1;lp<=l;lp++) 
		      new += mprob[j*L2+l*L+lp];
             pcost[i*L+l] = new;
           }
     free (mprob);

                /* find best factors. this takes O(wl^2) time in the worst 
		   case but the average should be closer to O(k*wl) */

     mbest = malloc ((wl+1) * (K+2) * sizeof (double));
     ibest = malloc ((wl+1) * (K+2) * sizeof (int));
     dinitial[0] = createMasks (K+1,m);
     for (k=1;k<=K;k++) dinitial[k] = dinitial[0] + k*maskSize(m);
     *wlen = 0; best = BK;
     for (l = L-1; l > 1; l--)
        { if (best < 1/(double)l) break;  /* cannot win */
          for (i=1;i<=wl;i++) mbest[i*(K+2)+(0)] = 0.0;
          for (k=1;k<=K+1;k++) 
             for (i=max(1,wl-k*l-(k-1)*tr+1);i<=wl;i++)
		 mbest[i*(K+2)+k] = 1.0;
          for (k=1; k<=K+1; k++)
             for (i=wl-k*l-(k-1)*tr; i>0; i--)
		 { cost = (pcost[i*L+l] < l+1) ?
		            pcost[i*L+l]/(l-pcost[i*L+l]+1) : 1.0;
		   if (cost > 1.0) cost = 1.0;
                   if (k==1) this = cost;
		   else this = 1-(1-cost)*(1-mbest[(i+l+tr)*(K+2)+(k-1)]);
                   ibest[i*(K+2)+k] = i;
                   if (mbest[(i+1)*(K+2)+k] < this)
                      { this = mbest[(i+1)*(K+2)+k];
                        ibest[i*(K+2)+k] = ibest[(i+1)*(K+2)+k];
                      }
                   mbest[i*(K+2)+k] = this;
                 }
          if (mbest[1*(K+2)+(K+1)] < best)  /* put the results in output */
             { *wlen = l;
               i = 1;
               for (k=K+1;k>=1;k--)
                  { COPY(dinitial[K+1-k],initials[ibest[i*(K+2)+k]],m);
                    i = ibest[i*(K+2)+k]+l+tr;
                  }
               best = mbest[1*(K+2)+(K+1)];
             }
        }
     free (pcost); free (mbest); free (ibest);
     free (initials[0]); free (initials);

     if (best >= BK) { *wlen = 0; return 1.0; }

	/* there is a promising k+1 splitting. select the appropriate states */

	/* compute dactive and dfinal */

     tmp = createMask(m);
     dfinal[0] = createMasks (K+1,m);
     dactive[0] = createMasks (K+1,m);
     for (k=0;k<=K;k++)
	{ dfinal[k] = dfinal[0] + k*maskSize(m);
	  dactive[k] = dactive[0] + k*maskSize(m);
	  COPY(dactive[k],dinitial[k],m);
	  for (i=0;i<*wlen;i++)
             { ZERO(tmp,m);
               for (j=0;j<m;j++)
                  if (ISSET(dactive[k],j)) OR (tmp,trans[j],m);
               OR(dactive[k],tmp,m);
	     }
	  COPY(dfinal[k],tmp,m);
        }
     free (tmp);
     return best;
   }

void eregularFreeScan (eregularScanData *scan)

   { regularFreeScan (scan->edata); 
     free (scan->V3);
     free (scan);
   }

                /* loads masks for fast bwd/fwd */

eregularScanData *eregularLoadFast (int m, Mask *trans, int wlen, 
				   Mask *B, Mask *active, Mask *initial, 
				   Mask *final, int k, int type, int **map)

        /* findBest ensures that this always fits in a single mask */

   { eregularScanData *scan = malloc (sizeof(eregularScanData));
     Mask uactive,uinitial,ufinal;
     int i,j,width;
     Mask *redtrans,redinitial,redfinal;
     Mask *revtrans,revinitial,revfinal;
     Mask redB[256];

     scan->k = k;
     scan->type = type;
     if (type != KPLUS1)
        { scan->edata = regularLoadFast(m,trans,wlen,B,*active,*initial,
					*final,map);
	}
     else
	{ uactive = ZERO(createMask(m),m);
	  uinitial = ZERO(createMask(m),m);
	  ufinal = ZERO(createMask(m),m);
	  for (i=0;i<=k;i++)
	     { OR(uactive,active[i],m);
	       OR(uinitial,initial[i],m);
	       OR(ufinal,final[i],m);
	     }
          scan->edata = regularLoadFast(m,trans,wlen,B,uactive,uinitial,ufinal,
					map);
	  free (uactive); free(uinitial); free(ufinal);
	}
     scan->dm = scan->edata->dm;
     scan->V3 = malloc (2*(k+1)*sizeof(mask));
     return scan;
   }

eregularData *eregularPreproc (byte *pat, Tree *tree, Tree **pos, int k)

   { eregularData *P = malloc (sizeof(eregularData));
     int i,c,wlen,wlen1,slices;
     Mask *trans,*rtrans;
     Mask *active,active1,*initial,initial1,*final,final1;
     Mask S[256],A;
     double best1,best2;
     int j,K,*beg,*end,*map;

	/* allocate and load the masks */

     P->m = 0;
     regularLength(tree,&P->m);
     P->initial = createMask(P->m);
     P->final = createMask(P->m);
     P->B[0] = createMasks (256,P->m);
     S[0] = createMasks (256,P->m);
     for (c=0; c<256; c++)
         { P->B[c] = ZERO(P->B[0]+c*maskSize(P->m),P->m);
           S[c] = ZERO(S[0]+c*maskSize(P->m),P->m);
         }
     A = ZERO(createMask(P->m),P->m);
     
     regularLoadMasks (pat,strlen(pat),P->m,tree,pos,P->B,S,A,
                       &trans,P->initial,P->final);

                /* make rtrans = reverse of trans, O(m^2/w + m^2) time */

     rtrans = malloc (P->m*sizeof(Mask));
     rtrans[0] = createMasks(P->m,P->m);
     for (i=0;i<P->m;i++)
        rtrans[i] = ZERO(rtrans[0]+i*maskSize(P->m),P->m);
     for (i=0;i<P->m;i++)
        for (j=0;j<P->m;j++)
           if (ISSET(trans[i],j)) SET(rtrans[j],i);

	/* determine a subset of the states for the search */
     active = malloc ((k+1)*sizeof(Mask*));
     initial = malloc ((k+1)*sizeof(Mask*));
     final = malloc ((k+1)*sizeof(Mask*));
     best1 = 1.3 * (k+1) * 
		     regularFindBest (tree,trans,rtrans,P->m,P->B,P->initial,
				      P->final,
				      &wlen1,&active1,&initial1,&final1,k);
     best2 = findBestMulti (tree,trans,rtrans,P->m,P->B,P->initial,P->final,
			    &wlen,active,initial,final,k);
     if ((best1 <= best2) || !wlen)  /* we prefer suffix */
	{ K = 0; wlen = wlen1; P->type = wlen ? BWD : FWD; 
	  active[0] = active1;
	  initial[0] = initial1;
	  final[0] = final1;
	}
     else
	{ P->type = KPLUS1; K = k;
          free (active1); free (initial1); free (final1);
	}

        /* fourth, create the verification data, common proc */

     P->k = k;
     slices = (P->m+OptDetWidth-1)/OptDetWidth;
     P->width = (P->m+slices-1)/slices;                                         
     regularMakeDet (P->width,trans,P->m,&P->dftransV);
     regularMakeDet (P->width,rtrans,P->m,&P->dbtransV);

	/* create searching subexpression according to the type of search */
     P->revMap = malloc(P->m*sizeof(int));
     switch (detClass(tree,active[0],P->m))
	{ case SIMPLE:
	     beg = malloc ((k+1)*sizeof(int));
	     P->match = ZEROS;
	     for (j=0;j<=k;j++)
	        { for (i=0;i<P->m;i++) if (ISSET(active[j],i)) break;
		  beg[j] = i;
		  if (wlen) P->match |= ONE<<i; 
	          for (i++;i<P->m;i++) if (!ISSET(active[j],i)) break;
		  if (!wlen) P->match |= ONE<<i;
		}
	     for (i=0;i<P->m;i++) P->revMap[i] = i;
	     P->dm = P->m;
	     P->scanData = esimpleLoadFast(wlen,k,P->type,P->B,beg);
	     free (beg);
	     P->scanText = (escanProc)esimpleScan;
	     P->scanFree = (freeScanProc)esimpleFreeScan;
	     break;
	  case EXTENDED:
	     beg = malloc ((k+1)*sizeof(int));
	     P->match = ZEROS;
	     for (j=0;j<=k;j++)
	        { for (i=0;i<P->m;i++) if (ISSET(active[j],i)) break;
		  beg[j] = i;
		  if (wlen) P->match |= ONE<<i; 
	          for (i++;i<P->m;i++) if (!ISSET(active[j],i)) break;
		  end[j] = i;
		  if (!wlen) P->match |= ONE<<i;
		}
	     for (i=0;i<P->m;i++) P->revMap[i] = i;
	     P->dm = P->m;
	     P->scanData = eextendedLoadFast(wlen,k,P->type,P->B,S,A,beg,end);
	     free (beg);
	     P->scanText = (escanProc)eextendedScan;
	     P->scanFree = (freeScanProc)eextendedFreeScan;
	     break;
	  case REGULAR:
	     P->scanData = eregularLoadFast(P->m,trans,wlen,P->B,active,
					    initial,final,k,P->type,&map);
	     for (i=0;i<P->m;i++) P->revMap[i] = -1;
             for (i=0;i<P->m;i++)
                 for (j=0;j<=K;j++)
                     if (ISSET(active[j],i)) P->revMap[map[i]] = i;
             free (map);                                                        
	     P->dm = ((eregularScanData*)P->scanData)->dm;
	     P->scanText = (escanProc)eregularScan;
	     P->scanFree = (freeScanProc)eregularFreeScan;
	     break;
        }

     free (S[0]); free (A);
     free (initial[0]); free (final[0]); free (active[0]);
     free (initial); free (final); free (active);
     free (trans[0]); free (trans);
     free (rtrans[0]); free (rtrans);
     P->V1 = createMasks (4,P->m);
     P->V2 = malloc (2*(k+1)*sizeof(Mask));
     P->V2[0] = createMasks (2*(k+1),P->m);
     for (i=0;i<2*(k+1);i++)
	P->V2[i] = P->V2[0] + i*maskSize(P->m);
     return P;
   }

void eregularFree (eregularData *P)

	/* Frees P */

   { P->scanFree (P->scanData);
     free (P->dftransV[0][0]); 
     free (P->dftransV[0]); free (P->dftransV);
     free (P->dbtransV[0][0]); 
     free (P->dbtransV[0]); free (P->dbtransV);
     free (P->initial); free (P->final);
     free (P->B[0]);
     free (P->revMap);
     free (P->V1); free (P->V2[0]); free (P->V2);
     free (P);
   }

static byte *fwdCheck (byte *ptr, byte *top, 
		       Mask **dtrans, Mask *B, Mask final, int m, 
		       int width, Mask *current, Mask *ocurrent, Mask T,
		       Mask tmp, Mask otmp, Mask ntmp, int *K)

		/* forward check.
		   receives maximal error level in *K and returns there
                   the best k found, being the return value the last
                   character considered in that case (NULL if nothing
                   can be found with *K errors)  */
		/* assumes k>0 and that ptr needs to be shifted before 
		   reading */

   { int i,f,j,n;
     int slices = (m+width-1)/width;
     int r,k = *K;
     byte *ret = NULL;
     Mask Boc,Bc;

           /* first remove trivial case of active final state */

     COPY(tmp,current[0],m);
     if (!ISZERO(AND(tmp,final,m),m))
        { *K = 0;
          while (*K <= k)
             { if (recCheckRightContext(ptr+1,top+1)) return ptr;
               if (ptr == top) return NULL; ptr++;
               if (OptIns) (*K)++; else return NULL;
             }
          return NULL;
        }

		/* fill initial masks */
     COPY(ocurrent[0],current[0],m);
     for (i=1;i<=k;i++)
	{ COPY(tmp,current[i-1],m);
	  if (OptDel)
	     { f = 0;
	       for (r=0;r<slices;r++)
	          { OR(tmp,dtrans[r][SLICE(current[i-1],f,width)],m);
	            f += width;
	            if ((f+width-1)/W != f/W) f += W - (f % W);
	          }
	       COPY(current[i],tmp,m);
               COPY(ocurrent[i],tmp,m);
			/* here it is possible that k >= m */
               if (!ISZERO(AND(tmp,final,m),m))
                  { if (recCheckRightContext(ptr+1,top+1))
                       { *K = i; k = i-1; ret = ptr; }
                  }
	      }
	}

           /* now process until automaton dies. each time it finds
              P with i<=k errors, annote i and pos and reduce k */
	
     n = 0;
     while (true)
	{ if (ptr == top) return ret;  /* end of buffer */
	  ptr++; n++;
	  Boc = Bc; Bc = B[*ptr]; /* new char */
	     /* now compute new state in ntmp */
	  ZERO (ntmp,m);
	  f = 0;
	  for (r=0;r<slices;r++)
	     { OR(ntmp,dtrans[r][SLICE(current[0],f,width)],m);
	       f += width;
	       if ((f+width-1)/W != f/W) f += W - (f % W);
	     }
    	  AND(ntmp,Bc,m);
          if (!ISZERO(AND(COPY(tmp,ntmp,m),final,m),m)) 
				/* found with zero errors, return */
             { if (recCheckRightContext(ptr+1,top+1))
                  { *K = 0; return ptr; }
	     }
	  COPY(otmp,current[0],m);
	  COPY(current[0],ntmp,m);
	  for (i=1;i<=k;i++)
	     { ZERO(tmp,m);
	       if (OptDel) 
		  { f = 0;
	            for (r=0;r<slices;r++)
	               { OR(tmp,dtrans[r][SLICE(ntmp,f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
		  }
	       COPY(ntmp,tmp,m);
	       if (OptIns) OR(ntmp,otmp,m);
	       if (OptSubs)
		  { f = 0;
	            for (r=0;r<slices;r++)
	               { OR(ntmp,dtrans[r][SLICE(otmp,f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
		  }
	       f = 0; ZERO(tmp,m);
	       for (r=0;r<slices;r++)
	          { OR(tmp,dtrans[r][SLICE(current[i],f,width)],m);
	            f += width;
	            if ((f+width-1)/W != f/W) f += W - (f % W);
	          }
	       AND(tmp,Bc,m);
               OR (ntmp,tmp,m);
	       if (OptTransp && (n >= 2))
		  { ZERO(T,m);
	            f = 0;
	            for (r=0;r<slices;r++)
	               { OR(T,dtrans[r][SLICE(ocurrent[i-1],f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
		    AND (T,Bc,m);
		    ZERO (tmp,m);
	            f = 0;
	            for (r=0;r<slices;r++)
	               { OR(tmp,dtrans[r][SLICE(T,f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
		    AND (tmp,Boc,m);
		    OR (ntmp,tmp,m);
		  }
	       COPY(tmp,ntmp,m);
               if (!ISZERO(AND(tmp,final,m),m)) /* match found, be stricter */
                  { if (recCheckRightContext(ptr+1,top+1))
                       { ret = ptr;
                         do { i--; COPY(tmp,current[i],m); }
                         while ((i >= 0) && !ISZERO(AND(tmp,final,m),m));
                         *K = i+1;
                         if (i == -1) return ptr; /* found with 0 errors */
                         k = i;
                       }
	          }
	       COPY(ocurrent[i-1],otmp,m);
	       COPY(otmp,current[i],m);
	       COPY(current[i],ntmp,m);
	     }
	  if (ISZERO(ntmp,m)) return ret;  /* automaton died, return best */
	}
   }

static byte *bwdCheck (byte *ptr, byte *top, 
		       Mask **dtrans, Mask *B, Mask final, int m, 
		       int width, Mask *current, Mask *ocurrent, Mask T,
		       Mask tmp, Mask otmp, Mask ntmp, int *K)

		/* backward check.
		   receives maximal error level in *K and returns there
                   the best k found, being the return value the last
                   character considered in that case (NULL if nothing
                   can be found with *K errors)  */
		/* assumes k>0 and that ptr needs to be shifted before 
		   reading */

   { int i,f,j,n;
     int slices = (m+width-1)/width;
     int r,k = *K;
     byte *ret = NULL;
     Mask Boc,Bc;

           /* first remove trivial case of active final state */

     COPY(tmp,current[0],m);
     if (!ISZERO(AND(tmp,final,m),m))
        { *K = 0;
          while (*K <= k)
             { if (recCheckLeftContext(ptr,top)) return ptr;
               if (ptr == top) return NULL; ptr--;
               if (OptIns) (*K)++; else return NULL;
             }
          return NULL;
        }

		/* fill initial masks */
     COPY(ocurrent[0],current[0],m);
     for (i=1;i<=k;i++)
	{ COPY(tmp,current[i-1],m);
	  if (OptDel)
	     { f = 0;
	       for (r=0;r<slices;r++)
	          { OR(tmp,dtrans[r][SLICE(current[i-1],f,width)],m);
	            f += width;
	            if ((f+width-1)/W != f/W) f += W - (f % W);
	          }
	       COPY(current[i],tmp,m);
               COPY(ocurrent[i],tmp,m);
			/* here it is possible that k >= m */
               if (!ISZERO(AND(tmp,final,m),m))
                  { if (recCheckLeftContext(ptr,top))
                       { *K = i; k = i-1; ret = ptr; }
                  }
	      }
	}

           /* now process until automaton dies. each time it finds
              P with i<=k errors, annote i and pos and reduce k */
	
     n = 0;
     while (true)
	{ if (ptr == top) return ret;  /* end of buffer */
	  ptr--; n++;
	  Boc = Bc; Bc = B[*ptr]; /* new char */
	     /* now compute new state in ntmp */
	  ZERO (ntmp,m);
	  f = 0;
	  COPY(tmp,current[0],m);
	  AND(tmp,Bc,m);
	  for (r=0;r<slices;r++)
	     { OR(ntmp,dtrans[r][SLICE(tmp,f,width)],m);
	       f += width;
	       if ((f+width-1)/W != f/W) f += W - (f % W);
	     }
          if (!ISZERO(AND(COPY(tmp,ntmp,m),final,m),m)) 
				/* found with zero errors, return */
             { if (recCheckLeftContext(ptr,top))
                  { *K = 0; return ptr; }
	     }
	  COPY(otmp,current[0],m);
	  COPY(current[0],ntmp,m);
	  for (i=1;i<=k;i++)
	     { ZERO(tmp,m);
	       if (OptDel) 
		  { f = 0;
	            for (r=0;r<slices;r++)
	               { OR(tmp,dtrans[r][SLICE(ntmp,f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
		  }
	       COPY(ntmp,tmp,m);
	       if (OptIns) OR(ntmp,otmp,m);
	       if (OptSubs)
		  { f = 0;
	            for (r=0;r<slices;r++)
	               { OR(ntmp,dtrans[r][SLICE(otmp,f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
		  }
	       f = 0;
	       COPY(tmp,current[i],m);
	       AND(tmp,Bc,m);
	       for (r=0;r<slices;r++)
	          { OR(ntmp,dtrans[r][SLICE(tmp,f,width)],m);
	            f += width;
	            if ((f+width-1)/W != f/W) f += W - (f % W);
	          }
	       if (OptTransp && (n >= 2))
		  { ZERO(T,m);
	            f = 0;
	            COPY(tmp,ocurrent[i-1],m);
	            AND(tmp,Bc,m);
	            for (r=0;r<slices;r++)
	               { OR(T,dtrans[r][SLICE(tmp,f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
	            f = 0;
		    AND(T,Boc,m);
	            for (r=0;r<slices;r++)
	               { OR(ntmp,dtrans[r][SLICE(T,f,width)],m);
	                 f += width;
	                 if ((f+width-1)/W != f/W) f += W - (f % W);
	               }
		  }
	       COPY(tmp,ntmp,m);
               if (!ISZERO(AND(tmp,final,m),m)) /* match found, be stricter */
                  { if (recCheckLeftContext(ptr,top)) 
                       { ret = ptr;
                         do { i--; COPY(tmp,current[i],m); }
                         while ((i >= 0) && !ISZERO(AND(tmp,final,m),m));
                         *K = i+1;
                         if (i == -1) return ptr; /* found with 0 errors */
                         k = i;
                       }
	          }
	       COPY(ocurrent[i-1],otmp,m);
	       COPY(otmp,current[i],m);
	       COPY(current[i],ntmp,m);
	     }
	  if (ISZERO(ntmp,m)) return ret;  /* automaton died, return best */
	}
   }

static bool checkMatch1 (eregularData *P, byte *pos, byte *beg, byte *end)

        /* checkMatch that does not care about a transposition where the
           pattern is divided */

   { byte *ptrb,*ptrf;          /* ptr to match of selected subpattern */
     int i,kb,kf;

        /* check which states have matched and test one by one */
     for (i=0;i<P->dm;i++)
         if ((P->match & (ONE<<i)) && (P->revMap[i] != -1))
            if (P->type == FWD)
	       {   /* first part: backward check */
                 kb = P->k;
		 ZERO(*P->V2,P->m);
		 SET(*P->V2,P->revMap[i]);
                 ptrb = bwdCheck (pos,beg,P->dbtransV,P->B,P->initial,P->m,
		      P->width,P->V2,P->V2+(P->k+1),P->V1,P->V1+maskSize(P->m),
		      P->V1+2*maskSize(P->m),P->V1+3*maskSize(P->m),&kb);
                 if (ptrb == NULL) continue;
	           /* second part: forward check */
                 kf = P->k-kb;
		 ZERO(*P->V2,P->m);
		 SET(*P->V2,P->revMap[i]);
                 ptrf = fwdCheck (pos-1,end-1,P->dftransV,P->B,P->final,P->m,
		      P->width,P->V2,P->V2+(P->k+1),P->V1,P->V1+maskSize(P->m),
		      P->V1+2*maskSize(P->m),P->V1+3*maskSize(P->m),&kf);
     	         if (ptrf == NULL) return false;
                 return true;   /* match in [ptrb,ptrf] */
	       }
            else
	       {   /* first part: forward check */
                 kf = P->k;
		 ZERO(*P->V2,P->m);
		 SET(*P->V2,P->revMap[i]);
                 ptrf = fwdCheck (pos,end-1,P->dftransV,P->B,P->final,P->m,
		      P->width,P->V2,P->V2+(P->k+1),P->V1,P->V1+maskSize(P->m),
		      P->V1+2*maskSize(P->m),P->V1+3*maskSize(P->m),&kf);
                 if (ptrf == NULL) continue;
	           /* second part: backward check */
                 kb = P->k-kf;
		 ZERO(*P->V2,P->m);
		 SET(*P->V2,P->revMap[i]);
                 ptrb = bwdCheck (pos+1,beg,P->dbtransV,P->B,P->initial,P->m,
		      P->width,P->V2,P->V2+(P->k+1),P->V1,P->V1+maskSize(P->m),
		      P->V1+2*maskSize(P->m),P->V1+3*maskSize(P->m),&kb);
     	         if (ptrb == NULL) return false;
                 return true;   /* match in [ptrb,ptrf] */
	       }
     return false;
   }

static bool checkMatch (eregularData *P, int dummy, 
			byte *pos, byte **beg, byte **end)

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

     ret = checkMatch1(P,pos,rbeg,rend);
     if (ret) { *beg = obeg; *end = oend; return true; }

       /* if the pattern is really split and transpositions are permitted
          a delicate problem may appear if the solution needs to transpose
          the characters at the split point. we solve the problem in a
          rather brutal but effective way */

     if (OptTransp && (pos < rend) && (pos-1 >= rbeg))
        { byte c = *pos;
          *pos = pos[-1];
          pos[-1] = c;
          ret = checkMatch1(P,pos,rbeg,rend);
          pos[-1] = *pos;
          *pos = c;
          if (ret) { *beg = obeg; *end = oend; return true; }
        }

     return false;
   }

static bool fwdScan11 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B,
		       mask initial, register mask final)
 
   { register byte *pos,*top;	/* current and final text pointers */
     register int c;
     register mask current0,current1;
     register mask ocurrent0,otmp;
     register mask *dtrans0 = dtrans[0];
     register mask Bc,Boc;

     pos = *beg; top = *end;
     while (pos < top)
	{ if ((c = *pos++) == OptRecChar) continue;
		/* fill initial masks */
          current0 = ocurrent0 = initial;
          current1 = current0 | dtrans0[current0];
             /* now traverse the text */
	  Bc = B[c];  /* take new char */
	     /* compute new state in tmp */
	  current1 = current0 | dtrans0[current0] | (dtrans0[current1] & Bc);
	  current0 = dtrans0[current0] & Bc;
	  current1 |= dtrans0[current0];
          Boc = Bc;
          while (pos < top)
	     { if (current1 & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = current1 & final;
	            if (checkMatch (P,0,pos,beg,end)) return true;
                  }
	       if ((c = *pos++) == OptRecChar) break;
	       Bc = B[c];  /* take new char */
	          /* now compute new state in tmp */
	       otmp = current0;
	       current0 = dtrans0[current0] & Bc;
	       current1 = otmp|dtrans0[current0|otmp]|(dtrans0[current1]&Bc)|
	                  (dtrans0[dtrans0[ocurrent0]&Bc]&Boc);
	       ocurrent0 = otmp;
               Boc = Bc;
	     }
	}
     return false;
   }

static bool fwdScan12 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B, 
		       mask initial, register mask final)
 
   { register byte *pos,*top;	/* current and final text pointers */
     register int c;
     register mask current0,current1,current2;
     register mask ocurrent0,ocurrent1,otmp,tmp;
     register mask *dtrans0 = dtrans[0];
     register mask Bc,Boc;

     pos = *beg; top = *end;
     while (pos < top)
	{ if ((c = *pos++) == OptRecChar) continue;
		/* fill initial masks */
          current0 = ocurrent0 = initial;
          current1 = ocurrent1 = current0 | dtrans0[current0];
          current2 = current1 | dtrans0[current1];
             /* now traverse the text */
	  Bc = B[c];  /* take new char */
	     /* compute new state in tmp */
	  current2 = current1 | dtrans0[current1] | (dtrans0[current2]&Bc);
	  current1 = current0 | dtrans0[current0] | (dtrans0[current1]&Bc);
	  current0 = dtrans0[current0] & Bc;
	  current1 |= dtrans0[current0];
	  current2 |= dtrans0[current1];
          Boc = Bc;
          while (pos < top)
	     { if (current2 & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = current2 & final;
	            if (checkMatch (P,0,pos,beg,end)) return true;
                  }
	       if ((c = *pos++) == OptRecChar) break;
	       Bc = B[c];  /* take new char */
	          /* now compute new state in tmp */
	       otmp = current0;
	       current0 = dtrans0[current0] & Bc;
	       tmp = otmp | dtrans0[current0|otmp] | (dtrans0[current1]&Bc) |
	             (dtrans0[dtrans0[ocurrent0]&Bc]&Boc);
	       ocurrent0 = otmp;
	       otmp = current1;
	       current1 = tmp;
	       current2 = otmp|dtrans0[current1|otmp]|(dtrans0[current2]&Bc)|
	                  (dtrans0[dtrans0[ocurrent1]&Bc]&Boc);
	       ocurrent1 = otmp;
               Boc = Bc;
	     }
	}
     return false;
   }

static bool fwdScan21 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B, 
		       mask initial, register mask final, register int width)
 
   { register mask aux;
     register byte *pos,*top;	/* current and final text pointers */
     register int c;
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask current0,current1;
     register mask ocurrent0,otmp;
     register mask *dtrans0 = dtrans[0];
     register mask *dtrans1 = dtrans[1];
     register mask Bc,Boc;

     pos = *beg; top = *end;
     while (pos < top)
	{ if ((c = *pos++) == OptRecChar) continue;
		/* fill initial masks */
          current0 = ocurrent0 = initial;
          current1 = current0 | dtrans0[current0&this] | dtrans1[current0>>width];
             /* now traverse the text */
	  Bc = B[c];  /* take new char */
	     /* compute new state */
	  aux = current0;
	  current0 = (dtrans0[current0&this] | dtrans1[current0>>width]) & Bc;
	  aux |= current0;
	  current1 = current0 | dtrans0[aux&this] | dtrans1[aux>>width] |
	             ((dtrans0[current1&this]|dtrans1[current1>>width])&Bc);
          Boc = Bc;
          while (pos < top)
	     { if (current1 & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = current1 & final;
	            if (checkMatch (P,0,pos,beg,end)) return true;
                  }
	       if ((c = *pos++) == OptRecChar) break;
	       Bc = B[c];  /* take new char */
	          /* now compute new state in tmp */
	       otmp = current0;
	       current0 = (dtrans0[current0&this] | dtrans1[current0>>width]) 
			  & Bc;
	       aux = otmp | current0;
	       current1 = otmp | dtrans0[aux&this]|dtrans1[aux>>width] |
	               ((dtrans0[current1&this]|dtrans1[current1>>width])&Bc);
	       aux = (dtrans0[ocurrent0&this]|dtrans1[ocurrent0>>width]) & Bc;
	       current1 |= (dtrans0[aux&this]|dtrans1[aux>>width]) & Boc;
	       ocurrent0 = otmp;
               Boc = Bc;
	     }
	}
     return false;
   }

static bool fwdScan22 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B, 
		       mask initial, register mask final, register int width)
 
   { register mask aux;
     register byte *pos,*top;	/* current and final text pointers */
     register int c;
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask current0,current1,current2;
     register mask ocurrent0,ocurrent1,otmp,tmp;
     register mask *dtrans0 = dtrans[0];
     register mask *dtrans1 = dtrans[1];
     register mask Bc,Boc;

     pos = *beg; top = *end;
     while (pos < top)
	{ if ((c = *pos++) == OptRecChar) continue;
		/* fill initial masks */
          current0 = ocurrent0 = initial;
          current1 = ocurrent1 = current0 |
	             dtrans0[current0&this] | dtrans1[current0>>width];
          current2 = current1 |
	             dtrans0[current1&this] | dtrans1[current1>>width];
           /* now traverse the text */
	  Bc = B[c];  /* take new char */
	     /* compute new state in tmp */
	  otmp = current0;
	  current0 = tmp = (dtrans0[current0&this] | dtrans1[current0>>width])
			   & Bc;
	  aux = tmp|otmp;
	  tmp = otmp | dtrans0[aux&this] | dtrans1[aux>>width];
	  tmp |= (dtrans0[current1&this] | dtrans1[current1>>width]) & Bc;
	  otmp = current1;
	  current1 = tmp;
	  aux = tmp|otmp;
	  current2 = otmp | dtrans0[aux&this] | dtrans1[aux>>width]
	            | ((dtrans0[current2&this]|dtrans1[current2>>width]) & Bc);
          Boc = Bc;
          while (pos < top)
	     { if (current2 & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = current2 & final;
	            if (checkMatch (P,0,pos,beg,end)) return true;
                  }
	       if ((c = *pos++) == OptRecChar) break;
	       Bc = B[c];  /* take new char */
	          /* now compute new state in tmp */
	       otmp = current0;
	       current0 = (dtrans0[current0&this] | dtrans1[current0>>width])
			  & Bc;
	       aux = current0|otmp;
	       tmp = otmp | dtrans0[aux&this] | dtrans1[aux>>width] |
	             ((dtrans0[current1&this]|dtrans1[current1>>width]) & Bc);
	       aux = (dtrans0[ocurrent0&this]|dtrans1[ocurrent0>>width]) & Bc;
	       tmp |= (dtrans0[aux&this]|dtrans1[aux>>width])&Boc;
	       ocurrent0 = otmp;
	       otmp = current1;
	       current1 = tmp;
	       aux = current1|otmp;
	       current2 = otmp | dtrans0[aux&this] | dtrans1[aux>>width] |
	             ((dtrans0[current2&this] | dtrans1[current2>>width]) & Bc);
	       aux = (dtrans0[ocurrent1&this]|dtrans1[ocurrent1>>width])&Bc;
	       current2 |= (dtrans0[aux&this]|dtrans1[aux>>width])&Boc;
	       ocurrent1 = otmp;
               Boc = Bc;
	     }
        }
     return false;
   }

static bool fwdScanrk (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, register mask **dtrans, 
		       register mask *B, mask initial, register mask final,
		       register int slices, register int width, register int k, 
		       register mask *current, register mask *ocurrent)
 
   { register int c,i,r;
     register mask tmp,otmp,aux,aux2,aux3,T;
     register byte *pos,*top;	/* current and final text pointers */
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask Bc,Boc;

     pos = *beg; top = *end;
     while (pos < top)
	{ if ((c = *pos++) == OptRecChar) continue;
		/* fill initial masks */
          tmp = current[0] = ocurrent[0] = initial;
          for (i=1;i<=k;i++)
	     { aux = tmp;
	       for (r=0;r<slices;r++)
	          { tmp |= dtrans[r][aux&this];
	            aux >>= width;
	          }
	       current[i] = ocurrent[i] = tmp;
	     }
             /* now traverse the text */
	  Bc = B[c];  /* take new char */
	     /* compute new state in tmp */
	  otmp = aux = current[0];
	  tmp = ZEROS;
	  for (r=0;r<slices;r++) 
	     { tmp |= dtrans[r][aux&this]; aux >>= width; }
	  tmp &= Bc;
	  current[0] = tmp;
	  for (i=1;i<=k;i++)
	     { aux = tmp|otmp; aux2 = current[i];
	       tmp = otmp;
	       for (r=0;r<slices;r++)
	          { tmp |= dtrans[r][aux&this] | (dtrans[r][aux2&this] & Bc);
	            aux >>= width; aux2 >>= width;
	          }
	       otmp = current[i];
	       current[i] = tmp;
             }
          Boc = Bc;
          while (pos < top)
	     { if (tmp & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = tmp & final;
	            if (checkMatch (P,0,pos,beg,end)) return true;
                  }
	       if ((c = *pos++) == OptRecChar) break;
	       Bc = B[c];  /* take new char */
	          /* now compute new state in tmp */
	       otmp = aux = current[0];
	       tmp = ZEROS;
	       for (r=0;r<slices;r++) 
		   { tmp |= dtrans[r][aux&this]; aux >>= width; }
	       tmp &= Bc;
	       current[0] = tmp;
	       for (i=1;i<=k;i++)
	          { aux = tmp|otmp; aux2 = current[i]; aux3 = ocurrent[i-1];
	            tmp = otmp; T = ZEROS;
	            for (r=0;r<slices;r++)
	               { tmp |= dtrans[r][aux&this]|(dtrans[r][aux2&this]&Bc);
	                 T |= dtrans[r][aux3&this];
	                 aux >>= width; aux2 >>= width; aux3 >>= width;
	               }
		    T &= Bc;
	            for (r=0;r<slices;r++) 
			{ tmp |= dtrans[r][T&this] & Boc; T >>= width;}
	            ocurrent[i-1] = otmp;
	            otmp = current[i];
	            current[i] = tmp;
                  }
               Boc = Bc;
	     }
	}
     return false;
   }

static bool kplus1Scan1 (byte **beg, byte **end, echeckProc checkMatch, 
		         eregularData *P, mask **dtrans, register mask *B, 
			 register mask *dtransin, register mask initial, 
			 register mask final, register int wlen, int k)
 
   { register mask icurrent,current;
     register byte *pos,*top;	/* current and final text pointers */
     register mask *dtrans0 = dtrans[0];
     register int j;
     int i;

     pos = *beg; top = *end;
     pos--;
     top -= wlen;
     while (pos < top)
        { if (!(current = dtransin[pos[wlen]]))  /* skip loop */
	     { pos += wlen; continue; }
          j = wlen-1;
	  while (true)
	     { icurrent = current & B[pos[j]];
	       current = dtrans0[icurrent];
	       if (!current) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
	            if (current & final)
		       { P->match = icurrent;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
	      }
	  pos += j;
        }
     return false;
   }

static bool kplus1Scan2 (byte **beg, byte **end, echeckProc checkMatch, 
		         eregularData *P, mask **dtrans, register mask *B, 
		         register mask *dtransin, register mask initial, 
			 register mask final, register int width, 
			 register int wlen, int k)
 
   { register mask icurrent,current;
     register byte *pos,*top;	/* current and final text pointers */
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask *dtrans0 = dtrans[0];
     register mask *dtrans1 = dtrans[1];
     register int j;
     int i;

     pos = *beg; top = *end;
     pos--;
     top -= wlen;
     while (pos < top)
        { if (!(current = dtransin[pos[wlen]]))  /* skip loop */
	     { pos += wlen; continue; }
	  j = wlen-1;
	  while (true)
	     {  /* compute transition */
	       icurrent = current & B[pos[j]];
	       current = dtrans0[icurrent&this]|dtrans1[icurrent>>width];
	       if (!current) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
	            if (current & final)
		       { P->match = icurrent;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
	      }
	  pos += j;
        }
     return false;
   }

static bool kplus1Scan (byte **beg, byte **end, echeckProc checkMatch, 
		        eregularData *P, mask **dtrans, register mask *B, 
			register mask *dtransin, register mask initial, 
			register mask final, register int slices, 
			register int width, register int wlen, int k)
 
   { register mask icurrent,current,tmp;
     register byte *pos,*top;	/* current and final text pointers */
     register int i,j;
     register mask this = (width == W) ? ONES : (ONE << width)-1;

     pos = *beg; top = *end;
     pos--;
     top -= wlen;
     while (pos < top)
        { if (!(current = dtransin[pos[wlen]]))  /* skip loop */
	     { pos += wlen; continue; }
	  j = wlen-1;
	  while (true)
	     {  /* compute transition */
	       icurrent = current & B[pos[j]]; 
	       tmp = icurrent; current = ZEROS;
	       for (i=0;i<slices;i++)
	           { current |= dtrans[i][tmp&this];
		     tmp >>= width;
	           }
	       if (!current) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
	            if (current & final)
		       { P->match = icurrent;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
	      }
	  pos += j;
        }
     return false;
   }

static bool bwdScan11 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B, 
		       register mask *dtransin, register mask initial, 
		       register mask final, register int wlen)
 
   { register byte *pos,*top;	/* current and final text pointers */
     register int j,c;
     register mask current0,current1;
     register mask ocurrent0,otmp,ocurrent1;
     register mask *dtrans0 = dtrans[0];
     register mask Bc,Boc;

		/* scan the text */
     pos = *beg; top = *end;
     pos--;
     wlen -= 1;
     top -= wlen;
     while (pos < top)
        {       /* fill initial masks (process one char now) */
	  c = pos[wlen]; /* take new character */
	  Boc = B[c]; /* take new character */
	  ocurrent0 = ocurrent1 = current1 = initial;
	  current0 = dtransin[c];
	  j = wlen-1;
	  while (true)
	     { Bc = B[pos[j]]; /* take new character */
		   /* compute transition */
	       otmp = current0;
	       current0 = dtrans0[current0 & Bc];
	       current1 = otmp|dtrans0[current0|otmp]|dtrans0[current1 & Bc]|
	                  dtrans0[dtrans0[ocurrent0 & Bc] & Boc];
	       ocurrent0 = otmp;
	       if (!current1) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
		    if (current1 & final)
		       { P->match = ocurrent1;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
               Boc = Bc;
	       ocurrent1 = current1;
	     }
	  pos += j;
        }
     return false;
   }

static bool bwdScan21 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B, 
		       register mask *dtransin, register mask initial, 
		       register mask final, register int width, 
		       register int wlen)
 
   { register mask aux,otmp;
     register byte *pos,*top;	/* current and final text pointers */
     register int j,c;
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask current0,current1;
     register mask ocurrent0,ocurrent1;
     register mask *dtrans0 = dtrans[0];
     register mask *dtrans1 = dtrans[1];
     register mask Bc,Boc;

		/* scan the text */
     pos = *beg; top = *end;
     pos--;
     wlen -= 2;
     top -= wlen;
     while (pos < top)
        {       /* fill initial masks (process one char now) */
	  c = pos[wlen]; /* take new character */
	  Boc = B[c]; /* take new character */
	  ocurrent0 = ocurrent1 = current1 = initial;
	  current0 = dtransin[c];
	  j = wlen-1;
	  while (true)
	     { Bc = B[pos[j]]; /* take new character */
		   /* compute transition */
	       otmp = current0;
	       current0 &= Bc;
	       current0 = dtrans0[current0&this] | dtrans1[current0>>width];
	       aux = current0|otmp;
	       current1 &= Bc;
	       current1 = otmp | dtrans0[aux&this] | dtrans1[aux>>width] |
	                  dtrans0[current1&this] | dtrans1[current1>>width];
	       ocurrent0 &= Bc;
	       aux = (dtrans0[ocurrent0&this]|dtrans1[ocurrent0>>width]) & Boc;
	       current1 |= dtrans0[aux&this] | dtrans1[aux>>width];
	       ocurrent0 = otmp;
	       if (!current1) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
		    if (current1 & final)
		       { P->match = ocurrent1;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
               Boc = Bc;
	       ocurrent1 = current1;
	     }
	  pos += j;
        }
     return false;
   }

static bool bwdScan12 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B, 
		       register mask *dtransin, register mask initial, 
		       register mask final, register int wlen)
 
   { register byte *pos,*top;	/* current and final text pointers */
     register int j,c;
     register mask current0,current1,current2;
     register mask ocurrent0,ocurrent1,ocurrent2,otmp,tmp;
     register mask *dtrans0 = dtrans[0];
     register mask Bc,Boc;

		/* scan the text */
     pos = *beg; top = *end;
     pos--;
     wlen -= 1;
     top -= wlen;
     while (pos < top)
        {       /* fill initial masks (process one char now) */
	  c = pos[wlen]; /* take new character */
	  Boc = B[c]; 
	  ocurrent0 = ocurrent1 = ocurrent2 = current1 = current2 = initial;
	  current0 = dtransin[c];
	  j = wlen-1;
	  while (true)
	     { Bc = B[pos[j]]; /* take new character */
		   /* compute transition */
	       otmp = current0;
	       current0 = dtrans0[current0 & Bc];
	       tmp = otmp | dtrans0[current0|otmp] | dtrans0[current1 & Bc] |
	             dtrans0[dtrans0[ocurrent0 & Bc] & Boc];
	       ocurrent0 = otmp;
	       otmp = current1;
	       current1 = tmp;
	       current2 = otmp|dtrans0[current1|otmp]|dtrans0[current2 & Bc] |
	                  dtrans0[dtrans0[ocurrent1 & Bc] & Boc];
	       ocurrent1 = otmp;
	       if (!current2) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
		    if (current2 & final)
		       { P->match = ocurrent2;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
               Boc = Bc;
	       ocurrent2 = current2;
	     }
	  pos += j;
        }
     return false;
   }

static bool bwdScan22 (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, mask **dtrans, register mask *B, 
		       register mask *dtransin, register mask initial, 
		       register mask final, register int width, 
		       register int wlen)
 
   { register mask aux,otmp,tmp,tmp2;
     register byte *pos,*top;	/* current and final text pointers */
     register int j,c;
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask current0,current1,current2;
     register mask ocurrent0,ocurrent1,ocurrent2;
     register mask *dtrans0 = dtrans[0];
     register mask *dtrans1 = dtrans[1];
     register mask Bc,Boc;

		/* scan the text */
     pos = *beg; top = *end;
     pos--;
     wlen -= 2;
     top -= wlen;
     while (pos < top)
        {       /* fill initial masks (process one char now) */
	  c = pos[wlen]; /* take new character */
	  Boc = B[c];
	  ocurrent0 = ocurrent1 = ocurrent2 = current1 = current2 = initial;
	  current0 = dtransin[c];
	  j = wlen-1;
	  while (true)
	     { Bc = B[pos[j]]; /* take new character */
		   /* compute transition */
	       otmp = current0;
	       current0 &= Bc;
	       current0 = dtrans0[current0&this] | dtrans1[current0>>width];
	       aux = current0|otmp;
	       tmp2 = current1 & Bc;
	       tmp = otmp | dtrans0[aux&this] | dtrans1[aux>>width] |
	             dtrans0[tmp2&this] | dtrans1[tmp2>>width];
	       tmp2 = ocurrent0 & Bc;
	       aux = (dtrans0[tmp2&this]|dtrans1[tmp2>>width]) & Boc;
	       tmp |= dtrans0[aux&this] | dtrans1[aux>>width];
	       ocurrent0 = otmp;
	       otmp = current1;
	       current1 = tmp;
	       aux = current1|otmp;
	       current2 &= Bc;
	       current2 = otmp | dtrans0[aux&this] | dtrans1[aux>>width] |
	             dtrans0[current2&this] | dtrans1[current2>>width];
	       aux = ocurrent1 & Bc;
	       aux = (dtrans0[aux&this] | dtrans1[aux>>width]) & Boc;
	       current2 |= dtrans0[aux&this] | dtrans1[aux>>width];
	       ocurrent1 = otmp;
	       if (!current2) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
	            if (current2 & final) 
		       { P->match = ocurrent2;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
               Boc = Bc;
	       ocurrent2 = current2;
	     }
	  pos += j;
        }
     return false;
   }

static bool bwdScanrk (byte **beg, byte **end, echeckProc checkMatch, 
		       eregularData *P, register mask **dtrans, 
		       register mask *B, register mask *dtransin,
		       register mask initial, register mask final, 
		       register int slices, register int width, 
		       register int wlen, register int k, 
		       register mask *current, register mask *ocurrent)
 
   { register mask tmp,otmp,aux,aux2,aux3,T;
     register byte *pos,*top;	/* current and final text pointers */
     register int i,j,r,c;
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask Boc,Bc;

		/* scan the text */
     pos = *beg; top = *end;
     pos--;
     wlen -= k;
     top -= wlen;
     while (pos < top)
        {       /* fill initial masks (process one char now) */
	  c = pos[wlen]; /* take new character */
	  Boc = B[c];
	  ocurrent[0] = initial;
	  current[0] = dtransin[c];
          for (i=1;i<=k;i++) current[i] = ocurrent[i] = initial;
	  j = wlen-1;
	  while (true)
	     { Bc = B[pos[j]]; /* take new character */
		   /* compute transition */
	       otmp = current[0];
	       aux = otmp & Bc;
	       tmp = ZEROS;
	       for (r=0;r<slices;r++)
	          { tmp |= dtrans[r][aux&this];
	            aux >>= width;
	          }
	       current[0] = tmp;
	       for (i=1;i<=k;i++)
	          { aux = tmp|otmp;
		    aux2 = current[i] & Bc;
		    aux3 = ocurrent[i-1] & Bc;
	            tmp = otmp; T = ZEROS;
	            for (r=0;r<slices;r++)
	                { tmp |= dtrans[r][aux&this] | dtrans[r][aux2&this];
	                  T |= dtrans[r][aux3&this];
		          aux >>= width; aux2 >>= width; aux3 >>= width;
	                }
		    T &= Boc;
	            for (r=0;r<slices;r++)
	               { tmp |= dtrans[r][T&this];
	                 T >>= width;
	               }
	            ocurrent[i-1] = otmp;
	            otmp = current[i];
	            current[i] = tmp;
                  }
	       if (!tmp) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
	            if (tmp & final) 
		       { P->match = otmp;
			 if (checkMatch (P,0,pos+1,beg,end)) return true;
		       }
	 		  /* else shift in 1 */
	            j = 1;
		    break;
	          }
               Boc = Bc;
	     }
	  pos += j;
        }
     return false;
   }

bool eregularScan (byte **beg, byte **end, echeckProc checkMatch, void *P,
		  eregularScanData *scan)
 
        /* Scans from *beg to *end-1 for P. Returns if it could find
           it. In case of returning true, *beg and *end are set to limit
           the first occurrence of P.
         */                                                                    

   { switch (scan->type)
	{ case KPLUS1:
	     switch (scan->edata->slices)
	        { case 1:
                  return kplus1Scan1 (beg,end,checkMatch,P,scan->edata->dtrans,
			            scan->edata->B,scan->edata->dtransin,
				    scan->edata->dinitial,scan->edata->dfinal,
			            scan->edata->wlen,scan->k);
	          case 2:
                  return kplus1Scan2 (beg,end,checkMatch,P,scan->edata->dtrans,
			            scan->edata->B,scan->edata->dtransin,
			            scan->edata->dinitial,scan->edata->dfinal,
			            scan->edata->width,scan->edata->wlen,
				    scan->k);
		  default:
                  return kplus1Scan (beg,end,checkMatch,P,scan->edata->dtrans,
			            scan->edata->B,scan->edata->dtransin,
			            scan->edata->dinitial,scan->edata->dfinal,
			            scan->edata->slices,scan->edata->width,
				    scan->edata->wlen,scan->k);
	        }
	  case FWD:
	     switch (scan->edata->slices)
	        { case 1:
		  switch (scan->k)
		     { case 1:
                       return fwdScan11 
			         (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dinitial,
				  scan->edata->dfinal);
		       case 2:
                       return fwdScan12 
			         (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dinitial,
				  scan->edata->dfinal);
		       default:
                       return fwdScanrk 
			      (beg,end,checkMatch,P,scan->edata->dtrans,
			       scan->edata->B,scan->edata->dinitial,
			       scan->edata->dfinal,scan->edata->slices,
			       scan->edata->width,scan->k,
		               scan->V3,scan->V3+(scan->k+1));
		     }
	          case 2:
		  switch (scan->k)
		     { case 1:
                       return fwdScan21 
				 (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dinitial,
				  scan->edata->dfinal,scan->edata->width);
		       case 2:
                       return fwdScan22 
				 (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dinitial,
				  scan->edata->dfinal,scan->edata->width);
		       default:
                       return fwdScanrk 
			      (beg,end,checkMatch,P,scan->edata->dtrans,
			       scan->edata->B,scan->edata->dinitial,
			       scan->edata->dfinal,scan->edata->slices,
			       scan->edata->width,scan->k,
		               scan->V3,scan->V3+(scan->k+1));
		     }
		  default:
                       return fwdScanrk 
			      (beg,end,checkMatch,P,scan->edata->dtrans,
			       scan->edata->B,scan->edata->dinitial,
			       scan->edata->dfinal,scan->edata->slices,
			       scan->edata->width,scan->k,
		               scan->V3,scan->V3+(scan->k+1));
	        }
	  case BWD: 
	     switch (scan->edata->slices)
	        { case 1:
		  switch (scan->k)
		     { case 1:
                       return bwdScan11 
			         (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dtransin,
				  scan->edata->dinitial,scan->edata->dfinal,
		                  scan->edata->wlen);
		       case 2:
                       return bwdScan12 
			         (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dtransin,
			          scan->edata->dinitial,scan->edata->dfinal,
		                  scan->edata->wlen);
		       default:
                       return bwdScanrk 
			      (beg,end,checkMatch,P,scan->edata->dtrans,
			       scan->edata->B,scan->edata->dtransin,
			       scan->edata->dinitial,scan->edata->dfinal,
		               scan->edata->slices,scan->edata->width,
			       scan->edata->wlen,
		               scan->k,scan->V3,scan->V3+(scan->k+1));
		     }
	          case 2:
		  switch (scan->k)
		     { case 1:
                       return bwdScan21 
				 (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dtransin,
			          scan->edata->dinitial,scan->edata->dfinal,
		                  scan->edata->width,scan->edata->wlen);
		       case 2:
                       return bwdScan22 
				 (beg,end,checkMatch,P,scan->edata->dtrans,
			          scan->edata->B,scan->edata->dtransin,
			          scan->edata->dinitial,scan->edata->dfinal,
		                  scan->edata->width,scan->edata->wlen);
		       default:
                       return bwdScanrk 
			      (beg,end,checkMatch,P,scan->edata->dtrans,
			       scan->edata->B,scan->edata->dtransin,
			       scan->edata->dinitial,scan->edata->dfinal,
		               scan->edata->slices,scan->edata->width,
			       scan->edata->wlen,
		               scan->k,scan->V3,scan->V3+(scan->k+1));
		     }
		  default:
                       return bwdScanrk 
			      (beg,end,checkMatch,P,scan->edata->dtrans,
			       scan->edata->B,scan->edata->dtransin,
			       scan->edata->dinitial,scan->edata->dfinal,
		               scan->edata->slices,scan->edata->width,
			       scan->edata->wlen,
		               scan->k,scan->V3,scan->V3+(scan->k+1));
	        }
	}
   }

bool eregularSearch (byte **beg, byte **end, eregularData *P)

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

   { return P->scanText (beg,end,(echeckProc)checkMatch,P,P->scanData);
   }
 
