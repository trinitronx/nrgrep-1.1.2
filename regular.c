
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

	/* Search a regular expression in a text */

#include "regular.h"
#include "simple.h"
#include "extended.h"

#define BF 0.65  /* threshold cost to prefer bwd over fwd */

double minCost (Tree *e, int m, int wl, int l, double *pcost, int *wlens,
		Mask mark, Mask markf, Mask arg, Mask iarg, Mask farg)

	/* Given e, pcost, wlens, mark*, finds the best necessary factor of
	   length l and delivers avg cost per window and the subautomaton 
	   selected in *arg. Assumes that *arg comes allocated and zeroed */

   { double c1,c2;
     int j,size;
     Mask tmp,tmpi,tmpf;
     switch (e->type)
       { case STR:
	    if (wlens[e->pos] < l) return l+1;
	    COPY(arg,mark+(e->pos*wl+l)*maskSize(m),m);
	    SET(ZERO(iarg,m),e->pos); 
	    COPY(farg,markf+(e->pos*wl+l)*maskSize(m),m);
	    return pcost[e->pos*wl+l];
	 case STAR: case QUESTION:
	    return l+1;
	 case PLUS:
	    return minCost(e->e1,m,wl,l,pcost,wlens,mark,markf,arg,iarg,farg);
	 case OOR:
	    tmp = ZERO(createMask(m),m);
	    tmpi = ZERO(createMask(m),m);
	    tmpf = ZERO(createMask(m),m);
	    c1 = minCost(e->e1,m,wl,l,pcost,wlens,mark,markf,arg,iarg,farg);
	    c2 = minCost(e->e2,m,wl,l,pcost,wlens,mark,markf,tmp,tmpi,tmpf);
	    OR(arg,tmp,m); OR(iarg,tmpi,m); COPY(farg,tmpf,m);
		/* this is a quick & dirty solution to the problem of 
		   limiting the number of bits to W. a better solution
		   would need dynprog over a third dimension: number of
		   bits used */
	    size = 0;
	    for (j=0;j<m;j++) if (ISSET(arg,j)) size++;
	    if (size > W) 
	       { ZERO(arg,m); ZERO(iarg,m); ZERO(farg,m); 
		 return l+1; 
	       }
	    return max(c1,c2);
	 case CONC:
	    tmp = ZERO(createMask(m),m);
	    tmpi = ZERO(createMask(m),m);
	    tmpf = ZERO(createMask(m),m);
	    c1 = minCost(e->e1,m,wl,l,pcost,wlens,mark,markf,arg,iarg,farg);
	    c2 = minCost(e->e2,m,wl,l,pcost,wlens,mark,markf,tmp,tmpi,tmpf);
	    if (c1 <= c2) return c1;
	    COPY(arg,tmp,m); free (tmp);
	    COPY(iarg,tmpi,m); free (tmpi);
	    COPY(farg,tmpf,m); free (tmpf);
	    return c2;
       }
   }

void regularComplete (int wlen, int m, Mask initial, Mask final, 
		      Mask *trans, Mask *rtrans,
		      Mask dinitial, Mask dactive, Mask dfinal,
		      Mask binitial, Mask bactive, Mask bfinal, 
		      Mask finitial, Mask factive, Mask ffinal, Mask **test)

		/* Given an automaton of size m, transitions trans,
		   initial,final states and reverse transitions rtrans, given
		   segment dinitial,dactive,dfinal of the autom; it writes the
		   segments binitial,bactive,bfinal for bwd verification and
		   segments finitial,factive,ffinal for fwd verification.
		   wlen>0 means bwd search of the d segment, otherwise is fwd
		*/

   { Mask noactive,tmp,tmp2;
     int i;
     COPY (binitial,initial,m); 
     COPY (ffinal,final,m);
     noactive = NOT(COPY(createMask(m),dactive,m),m); 
     *test = malloc (m*sizeof(Mask*));
     (*test)[0] = createMasks(m,m);
     for (i=0;i<m;i++)
	 { (*test)[i] = (*test[0]) + i*maskSize(m);
	   ZERO((*test)[i],m);
	 }
     if (wlen)
	{ COPY(finitial,dinitial,m); 
		/* bfinal = previous to dinitial. there can be some extra
		   states here */
	  ZERO(bfinal,m);
	  for (i=0;i<m;i++)
	     if (ISSET(dinitial,i)) 
		{ OR (bfinal,rtrans[i],m);
		  COPY ((*test)[i],rtrans[i],m);
		}
	}
     else
	{ COPY (bfinal,dfinal,m);
		/* finitial = next to dfinal. there can be some extra
		   states here */
	  ZERO(finitial,m);
	  for (i=0;i<m;i++)
	     if (ISSET(dfinal,i)) 
		{ OR (finitial,trans[i],m);
		  COPY((*test)[i],trans[i],m);
		}
	}
	/* compute active states by closure */
     if (ISZERO(finitial,m)) ZERO(factive,m); else COPY(factive,final,m);
     if (ISZERO(bfinal,m)) ZERO(bactive,m); else COPY(bactive,initial,m);
     tmp = createMask(m);
     tmp2 = createMask(m);
     if (wlen)
        { while (true)
	    { COPY(tmp,bactive,m);
	      for (i=0;i<m;i++)
	         if (ISSET(bactive,i))
	            OR(bactive,AND(COPY(tmp2,noactive,m),trans[i],m),m);
	      if (EQUAL(tmp,bactive,m)) break;
	    }
	  COPY(noactive,bactive,m);
	  NOT(noactive,m);
          while (true)
	    { COPY(tmp,factive,m);
	      for (i=0;i<m;i++)
	         if (ISSET(factive,i))
	            OR(factive,AND(COPY(tmp2,noactive,m),rtrans[i],m),m);
	      if (EQUAL(tmp,factive,m)) break;
	    }
	}
     else
        { while (true)
	    { COPY(tmp,factive,m);
	      for (i=0;i<m;i++)
	         if (ISSET(factive,i))
	            OR(factive,AND(COPY(tmp2,noactive,m),rtrans[i],m),m);
	      if (EQUAL(tmp,factive,m)) break;
	    }
	  COPY(noactive,factive,m);
	  NOT(noactive,m);
          while (true)
	    { COPY(tmp,bactive,m);
	      for (i=0;i<m;i++)
	         if (ISSET(bactive,i))
	            OR(bactive,AND(COPY(tmp2,noactive,m),trans[i],m),m);
	      if (EQUAL(tmp,bactive,m)) break;
	    }
	}
     free (tmp); free (tmp2); free (noactive);
	/* now purgue unwanted finitial,bfinal */
     AND (finitial,factive,m); AND (ffinal,factive,m);
     AND (bfinal,bactive,m); AND (bfinal,bactive,m);
     if (wlen)
	  for (i=0;i<m;i++) AND((*test)[i],bactive,m);
     else for (i=0;i<m;i++) AND((*test)[i],factive,m);
   }

static void regularGetMost (Mask *trans, int m,
			    Mask initial, Mask final,
			    Mask dactive, Mask dinitial, Mask dfinal)

	/* Gets all the states that fit in a computer word, trying to
	   minimize the probability of getting outside these states. This
	   is done in a very simple way: the prob. of reaching a state is
	   simply computed as inversely proportional to its shortest path
	   from the initial state */

   { Mask oreach,reach;
     int j;
		/* simple case: all fits */
     COPY (dinitial,initial,m);
     ZERO(dactive,m);
     if (m <= W)
	{ for (j=0;j<m;j++) SET(dactive,j);
	  COPY(dfinal,final,m);
	  return;
	}
		/* add states until they don't fit, O(w*m^2/w) = O(m^2) */
     reach = COPY(createMask(m),initial,m);
     oreach = createMask(m);
     while (ONEBITS(reach,m) <= W)
        { COPY(oreach,dactive,m);
	  NOT(oreach,m);
	  AND(oreach,reach,m);
	  OR(dactive,reach,m);
	  ZERO(reach,m);
	  for (j=0;j<m;j++)
	      if (ISSET(dactive,j)) OR (reach,trans[j],m);
	}
     COPY (dfinal,oreach,m);
     AND(OR(dfinal,final,m),dactive,m);
     free (oreach); free (reach);
   }

double regularFindBest (Tree *e, Mask *trans, Mask *rtrans, int m, 
			Mask *B, Mask initial, Mask final, int *wlen, 
			Mask *dactive, Mask *dinitial, Mask *dfinal, int K)

        /* Finds a best subexpression to search and puts in dactive,dinitial,
	   dfinal its states and in wlen the window length (0 if it recommends
	   forward scanning). Returns the avg number of chars inspected. K is 
	   a number to add to all the costs */

   { double *prob;    /* prob[i] = prob of arriving to state i from prev */
     double *pprob;   /* pprob[i,l] = prob matching a path of length l from
				     state i */
     mask *mark;      /* mark[i,l] = reachable nodes in l steps from i */
     mask *markf;     /* markf[i,l] = last reachable nodes in l steps from i */
     double *mprob;   /* mprob[i,l,l'] = prob of any factor of length l'
                         inside area mark[i,l] */
     double *pcost;    /* pcost[i,l] = avg cost per window in area mark[i,l] */
     int *wlens;      /* wlens[i] = shortest path from i to a final state */
     int wl,c,i,j,l,lp,wl2,bestl;
     Mask reach,reachi,reachf,tmp,tmpi,tmpf;
     double best,new;


		/* compute max wlen wl <= m, wlens[],
		   in O(wl*m^3/w) = O(m^4/w) time */

     reach = createMask(m);
     tmp = createMask(m);
     wlens = malloc (m * sizeof(int));
     for (i=0;i<m;i++)
	{ wlens[i] = 1; /* there is an extra initial state */
          SET(ZERO(reach,m),i);
          while (ISZERO(AND(COPY(tmp,reach,m),final,m),m))
	     { wlens[i]++;
	       ZERO(tmp,m);
	       for (j=0;j<m;j++)
	          if (ISSET(reach,j)) OR (tmp,trans[j],m);
	       OR(reach,tmp,m);
	     }
	}
     free (tmp); free (reach);
     wl = wlens[0]; if (wl > W) wl = W;
     wl2 = wl*wl; /* for simplicity */

	 	/* load letter probabilities in O(m*s) time */

     prob = malloc (m * sizeof (double));
     for (i=0;i<m;i++) 
         { prob[i] = 0.0;
           for (c=0;c<256;c++) 
	       if (ISSET(B[c],i)) prob[i] += letterProb[c];
	 }
     
		/* compute path probabilities in O(wl*m^2) time */
		/* mark reachable states, O(wl*m^3/w) time */
		/* pprob[i,l] includes already the prob. of reaching i
		   from its predecessor */

     pprob = malloc (m * wl * sizeof (double));
     mark = createMasks (m * wl, m);
     markf = createMasks (m * wl, m);
     for (i=0;i<m;i++) 
	{ pprob[i*wl+0] = 1.0;
	  ZERO(mark+(i*wl+0)*maskSize(m),m);
	  ZERO(markf+(i*wl+0)*maskSize(m),m);
	  pprob[i*wl+1] = prob[i];
	  SET(ZERO(mark+(i*wl+1)*maskSize(m),m),i);
	  SET(ZERO(markf+(i*wl+1)*maskSize(m),m),i);
	}
     for (l=2;l<wl;l++)
	for (i=0;i<m;i++)
	   { pprob[i*wl+l] = 0.0;
	     ZERO(markf+(i*wl+l)*maskSize(m),m);
	     for (j=0;j<m;j++)
		if (ISSET(trans[i],j)) 
		   { pprob[i*wl+l] += pprob[j*wl+(l-1)];
	  	     OR(markf+(i*wl+l)*maskSize(m),
			markf+(j*wl+(l-1))*maskSize(m),m);
		   }
	     pprob[i*wl+l] *= prob[i];
	     if (pprob[i*wl+l] > 1.0) pprob[i*wl+l] = 1.0;
	     COPY(mark+(i*wl+l)*maskSize(m),markf+(i*wl+l)*maskSize(m),m);
	     OR(mark+(i*wl+l)*maskSize(m),mark+(i*wl+(l-1))*maskSize(m),m);
	   }
     free (prob);
	      
		/* compute mprob in O(m^2*wl^2) time */

     mprob = malloc (m * wl * wl * sizeof (double));
     for (i=0;i<m;i++)
        for (l=0;l<wl;l++)
           mprob[i*wl2+l*wl+0] = 1.0;
     for (l=1;l<wl;l++)
        for (lp=1;lp<=l;lp++)
	   for (i=0;i<m;i++)
	      { mprob[i*wl2+l*wl+lp] = pprob[i*wl+lp];
		if (lp < l)
	           for (j=0;j<m;j++)
		      if (ISSET(trans[i],j))
			 mprob[i*wl2+l*wl+lp] = 1 -
			  (1-mprob[i*wl2+l*wl+lp])*(1-mprob[j*wl2+(l-1)*wl+lp]);
	      }
     free (pprob);

		/* compute pcost in O(wl^2*m) time */

     pcost = malloc (m * wl * sizeof (double));
     for (i=0;i<m;i++)
        for (l=0;l<wl;l++)
	   { new = K;
	     for (lp=0;lp<=l;lp++) new += mprob[i*wl2+l*wl+lp];
	     pcost[i*wl+l] = new;
	   }
     free (mprob);

		/* find the best choice, in O(wl*m^2/w) time */

     tmp = createMask(m); tmpi = createMask(m); tmpf = createMask(m);
     reach = createMask(m); reachi = createMask(m); reachf = createMask(m);
     best = BF; bestl = 0;
     for (l=wl-1;l>=1;l--)
	{ if (best*l <= 1.0) break;  /* no chance with this length */
	  ZERO(tmp,m); ZERO(tmpi,m); ZERO(tmpf,m);
	  new = minCost (e,m,wl,l,pcost,wlens,mark,markf,tmp,tmpi,tmpf);
	  if (new<l+1) 
	     { new = new/(l-new+1);
	       if (new < best)
	          { best = new; bestl = l;
	            COPY(reach,tmp,m); COPY(reachi,tmpi,m); COPY(reachf,tmpf,m);
	          }
	     }
	}
     free (tmp); free (tmpi); free (tmpf);
     free (mark); free (markf); free (pcost); free (wlens);

		/* now set the output variables */

     *dactive = reach; *dinitial = reachi; *dfinal = reachf;
     if (best < BF) 
	{ *wlen = bestl; 
	}
     else 
	{ *wlen = 0;
	  regularGetMost (trans,m,initial,final,*dactive,*dinitial,*dfinal);
	}
     return best;
   }

	/* compute number of bits necessary to hold the regexp */

void regularLength (Tree *e, int *L)

   { if (*L == 0) *L = 1;
     switch (e->type)
        { case STR:
             if (!e->eps) e->pos = (*L)++;
             break;
          case CONC: case OOR:
             regularLength(e->e1,L);
             regularLength(e->e2,L);
             break;
          case PLUS: case STAR: case QUESTION:
             regularLength(e->e1,L);
             break;
        }
   }

static void firstLast (Tree *e, int L) 

	/* compute firstPos and lastPos, initial and final positions */

  { switch (e->type) 
	{ case STR: 
	     e->firstPos = ZERO(createMask(L),L);
	     if (!e->eps) e->firstPos = SET(e->firstPos,e->pos);
             e->lastPos = COPY(createMask(L),e->firstPos,L);
             break;
          case STAR: case PLUS: case QUESTION:
	     firstLast (e->e1,L);
             e->firstPos = COPY(createMask(L),e->e1->firstPos,L);
             e->lastPos = COPY(createMask(L),e->e1->lastPos,L);
             break;
          case OOR: 
	     firstLast (e->e1,L); firstLast (e->e2,L);
             e->firstPos = OR(COPY(createMask(L),e->e1->firstPos,L),
			      e->e2->firstPos,L);
             e->lastPos = OR(COPY(createMask(L),e->e1->lastPos,L),
			     e->e2->lastPos,L);
             break;     
          case CONC:
	     firstLast (e->e1,L); firstLast (e->e2,L);
             if (e->e1->eps) 
		  e->firstPos = OR(COPY(createMask(L),e->e1->firstPos,L),
                                   e->e2->firstPos,L);
      	     else e->firstPos = COPY(createMask(L),e->e1->firstPos,L);
             if (e->e2->eps) 
		  e->lastPos = OR(COPY(createMask(L),e->e1->lastPos,L),
                                  e->e2->lastPos,L);
      	     else e->lastPos = COPY(createMask(L),e->e2->lastPos,L);
             break;   
        }
  }

static Mask follow (Tree *e, int pos, int L)

	/* computes the follow set for a given position pos, that is,
	   the set of states which can follow position pos */

  { Mask tmp1,tmp2;
    switch (e->type)
	{ case STR:
             return ZERO(createMask(L),L);
             break;
          case STAR: case PLUS:
             if (ISSET(e->e1->lastPos,pos)) 
		  return OR(follow(e->e1,pos,L),e->e1->firstPos,L);
             else return follow(e->e1,pos,L);
             break;
          case QUESTION: 
             return follow(e->e1,pos,L);
             break;
          case OOR: 
             if (ISSET(e->e1->maskPos,pos))  /* pos appears in left subtree */
		{ if (ISSET(e->e2->maskPos,pos))  /* also in right subtree */
		     { tmp1 = follow(e->e1,pos,L);
	               tmp2 = follow(e->e2,pos,L);
		       OR (tmp1,tmp2,L); free (tmp2);
		       return tmp1;
		     }
		  else /* pos does not appear in right subtree */
		       return follow(e->e1,pos,L);
      	        } 
	     else  /* pos does not appear in left subtree */
		{ if (ISSET(e->e2->maskPos,pos)) /* appears in right subtree */
	   	       return follow(e->e2,pos,L);
		  else  /* does not appear anywhere */
		       return ZERO(createMask(L),L);
      	        }
             break;
          case CONC: 
      	     tmp1 = createMask(L);
             if (ISSET(e->e1->lastPos,pos)) COPY(tmp1,e->e2->firstPos,L);
      	     else ZERO(tmp1,L);
             if (ISSET(e->e1->maskPos,pos))  /* pos appears in left subtree */
		{ if (ISSET(e->e2->maskPos,pos))  /* also in right subtree */
		     { tmp2 = follow(e->e1,pos,L);
		       OR (tmp1,tmp2,L); free (tmp2);
	               tmp2 = follow(e->e2,pos,L);
		       OR (tmp1,tmp2,L); free (tmp2);
		       return tmp1;
		     } 
		  else  /* pos does not appear in right subtree */
		     { tmp2 = follow(e->e1,pos,L);
		       OR (tmp1,tmp2,L); free (tmp2);
	  	       return tmp1;
		     }
                } 
	     else  /* pos does not appear in left subtree */
		{ if (ISSET(e->e2->maskPos,pos)) /* appears in right subtree */
		     { tmp2 = follow(e->e2,pos,L);
		       OR (tmp1,tmp2,L); free (tmp2);
		       return tmp1;
		     } 
		  else  /* does not appear anywhere */
		       return tmp1;
      	        }
             break;   
         }
    }

	/* parses and builds the automaton from regular expression str */


static void loadExtended (Tree *tree, Mask *B, Mask *S, Mask A)

                /* reads masks S,A from pattern where it behaves as an
		   extended pattern */

   { int c,l;
     switch (tree->type)
        { case STR:
            break;
          case PLUS: 
	    if (tree->e1->type == STR) /* add in S what we did to B */
               { l = tree->e1->pos;
                 for (c=0; c<256; c++)
                     if (ISSET(B[c],l)) SET(S[c],l);
	       }
            break;
          case QUESTION: 
	    if (tree->e1->type == STR) /* mark S */
               SET(A,tree->e1->pos);
            break;
          case STAR:  /* treat it as +? */
	    if (tree->e1->type == STR)
               { l = tree->e1->pos;
                 for (c=0; c<256; c++)
                     if (ISSET(B[c],l)) SET(S[c],l);
                 SET(A,l);
	       }
            break;
          case CONC: /* recurse */
            loadExtended (tree->e1,B,S,A);
            loadExtended (tree->e2,B,S,A);
            break;
        }
   }

	/* receives pat and its strlen m, the number of bits L to store the
	   regexp, trees e and pos[]. It fills the preallocated tables B,S,A,
	   and the automaton: transitions *trans (not preallocated), and
	   initial and final states (preallocated) */

void regularLoadMasks (byte *pat, int m, int L, Tree *e, Tree **pos,
		       Mask *B, Mask *S, Mask A,
                       Mask **trans, Mask initial, Mask final)

   { int i,j;

	/* compute the B,S,A tables */
     i = 0;
     while (i<m)
	 if (pos[i]) getAclass(pat,&i,B,pos[i]->pos);
	 else i++;
     if (OptRecChar != -1) ZERO(B[OptRecChar],L);
     loadExtended (e,B,S,A);

	/* compute maskPos */
     setMaskPos (e,L);
	/* compute firstPos and lastPos */
     firstLast (e,L);

	/* initial and final states */ 
     SET(ZERO(initial,L),0);
     COPY(final,e->lastPos,L);
    
	/* allocate the transitions */
     *trans = malloc (L*sizeof(Mask*));
     (*trans)[0] = createMasks (L,m);
     for (i=1;i<L;i++) 
         { (*trans)[i] = (*trans)[0] + i*maskSize(m);
           ZERO((*trans)[i],m);
	 }

	/* load the transitions. this is basically the follow set, but it is
	    "first" for the initial state */
     for (i=0;i<L;i++) 
	{ if (i==0) COPY((*trans)[i],e->firstPos,L); /* initial position */
          else  /* all other positions */
             { Mask tmp = follow(e,i,L);
	       COPY((*trans)[i],tmp,L);
	       free (tmp);
	     }
        }
   }

void regularReverseArrows (Mask *trans, int m, Mask initial, 
		  Mask final, Mask **rtrans, Mask *rinitial, Mask *rfinal)

	/* reverses all the arrows of the NFA, exchanges initial and
  	   final states. the eps closures must have been done already */

   { int i,j;
     *rinitial = COPY(createMask(m),final,m);
     *rfinal = COPY(createMask(m),initial,m);
     *rtrans = malloc (m*sizeof(Mask*));
     (*rtrans)[0] = createMasks(m,m);
     for (i=0;i<m;i++)
         { (*rtrans)[i] = (*rtrans)[0] + i*maskSize(m);
	   ZERO((*rtrans)[i],m);
	 }
     for (i=0;i<m;i++) 
        for (j=0;j<m;j++) 
	   if (ISSET(trans[i],j)) SET((*rtrans)[j],i);
   }

void regularRemapStates (Mask *trans, int m, Mask initial, Mask final,
		  Mask active, Mask *B, Mask **rtrans, Mask *rB, int *rm, 
		  int *width, int *slices, Mask *rinitial, Mask *rfinal,
		  int **Map)

	/* remaps the states so that only the active ones exist. it does
	   not check that the reduction is consistent. width is the width
	   that will be used for the deterministic tables (<= W) and the
	   mapping guarantees that no slice will be cut by a machine word.
	   A new initial state 0 is used (provided it is not already there).
	   map is a memory area of m integers where the mapping of states is
	   returned, if the caller wants to use it for something else.
	*/

   { int i,c,j,*map;
	
		/* first compute number of states and slices/width */
     *rm = 1;
     for (i=1;i<m;i++) 
	 if (ISSET(active,i)) (*rm)++;
     *slices = (*rm+OptDetWidth-1)/OptDetWidth;
     *width = (*rm+*slices-1)/(*slices);
		/* now map the states */
     *Map = map = malloc (m * sizeof(int));
     j = 0; map[0] = 0; /* valid be state 0 active or not */
     for (i=1;i<m;i++) 
	 if (ISSET(active,i)) 
	    { j++;
	      map[i] = j; 
	    }
		/* remap initial state, always zero */
     *rinitial = SET(ZERO(createMask(*rm),*rm),0);
		/* remap final states */
     *rfinal = ZERO(createMask(*rm),*rm);
     for (i=0;i<m;i++) 
         if (ISSET(final,i)) SET(*rfinal,map[i]);
		/* compute transitions. this also works if 0 is
		   already the initial state */
     *rtrans = malloc (*rm*sizeof(Mask*));
     (*rtrans)[0] = createMasks (*rm,*rm);
     for (i=0;i<*rm;i++) 
         { (*rtrans)[i] = (*rtrans)[0] + i*maskSize(*rm);
           ZERO((*rtrans)[i],*rm);
         }
     rB[0] = createMasks (256,*rm);
     for (c=1;c<256;c++) 
         { rB[c] = rB[0] + c*maskSize(*rm);
           ZERO(rB[c],*rm);
         }
		/* from the initial state */
     for (j=0;j<m;j++) 
         if (ISSET(initial,j)) SET((*rtrans)[0],map[j]);
		/* from the rest (can include initial state too!) */
     for (i=0;i<m;i++) 
        if (ISSET(active,i))
           for (j=0;j<m;j++) 
               if (ISSET(active,j))
	          if (ISSET(trans[i],j)) 
		     SET((*rtrans)[map[i]],map[j]);
		/* map characters in B */
     for (c=0;c<256;c++) ZERO(rB[c],*rm);
     for (j=0;j<m;j++) 
        if (ISSET(active,j))
          for (c=0;c<256;c++)
	     if (ISSET(B[c],j)) SET(rB[c],map[j]);
   }

void regularMakeDet (int width, Mask *trans, int m, Mask ***dtrans)

	/* builds a deterministic table for trans. it can handle multiwords
	   in the det table too, but this is compatible with monoword usage.
	   it receives the width of the deterministic table. the NFA is
	   built so that the slices are always inside a single word */

   { int slices;
     int i,j,f,w,b;
     mask m1,m2;
     int dm = 1 << width;
     slices = (m+width-1)/width;
	/* allocate and structure the memory */
     *dtrans = malloc (slices * sizeof(Mask**));
     (*dtrans)[0] = malloc (slices * dm * sizeof(Mask*));
     for (i=1;i<slices;i++) 
	 (*dtrans)[i] = (*dtrans)[0] + i*dm;
     (*dtrans)[0][0] = createMasks (slices*dm,m);
     for (i=0;i<slices;i++) 
	 { (*dtrans)[i][0] = (*dtrans)[0][0] + i*dm*maskSize(m);
	   for (j=1;j<dm;j++)
	       (*dtrans)[i][j] = (*dtrans)[i][0] + j*maskSize(m);
	 }
	/* fill with NFA */
     f = 0;
     for (i=0;i<slices;i++) 
         { w = (i < slices-1) ? width : m-f;
           ZERO((*dtrans)[i][0],m);
           for (b=0; b < w; b++)
               { if (b+1 < W)
		    { m1 = ~(1<<b); m2 = 1 << (b+1);
                      for (j=1<<b; j<m2; j++)
			 { COPY((*dtrans)[i][j],(*dtrans)[i][j&m1],m);
                           OR((*dtrans)[i][j],trans[f+b],m);
			 }
		    }
	         else
		    { m1 = ~(1<<b);
                      for (j=1<<b; j != 0; j++)
			  { COPY((*dtrans)[i][j],(*dtrans)[i][j&m1],m);
                            OR((*dtrans)[i][j],trans[f+b],m);
			  }
		    }
               }
	   f += width;
	 }
   }

void regularMakeDet1 (int width, Mask *trans, int m, mask ***dtrans)

	/* same as regularMakeDet but the masks are simple */

   { int slices;
     int i,j,f,w,b;
     mask m1,m2;
     int dm = 1 << width;
     slices = (m+width-1)/width;
	/* allocate and structure the memory */
     *dtrans = malloc (slices * sizeof(mask**));
     (*dtrans)[0] = malloc (slices * dm * sizeof(mask));
     for (i=1;i<slices;i++) 
	 (*dtrans)[i] = (*dtrans)[0] + i*dm;
	/* fill with NFA */
     f = 0;
     for (i=0;i<slices;i++) 
         { w = (i < slices-1) ? width : m-f;
           (*dtrans)[i][0] = ZEROS;
           for (b=0; b < w; b++)
               { if (b+1 < W)
		    { m1 = ~(1<<b); m2 = 1 << (b+1);
                      for (j=1<<b; j<m2; j++)
			  (*dtrans)[i][j] = (*dtrans)[i][j&m1] | *(trans[f+b]);
		    }
	         else
		    { m1 = ~(1<<b);
                      for (j=1<<b; j != 0; j++)
			  (*dtrans)[i][j] = (*dtrans)[i][j&m1] | *(trans[f+b]);
		    }
               }
	   f += width;
	 }
   }

void regularFreeScan (regularScanData *scan)

   { free (scan->dtrans[0]); 
     free (scan->dtrans);
     free (scan);
   }

                /* loads masks for fast bwd/fwd */

regularScanData *regularLoadFast (int m, Mask *trans, int wlen, Mask *B,
				  Mask active, Mask initial, Mask final, 
				  int **map)

        /* findBest ensures that this always fits in a single mask */

   { int i,c;
     regularScanData *scan = malloc (sizeof(regularScanData));
     Mask *redtrans,redinitial,redfinal;
     Mask *revtrans,revinitial,revfinal;
     mask tmp,this;
     Mask redB[256];

     scan->wlen = wlen;
     regularRemapStates (trans,m,initial,final,active,B,&redtrans,redB,
	   &scan->dm,&scan->width,&scan->slices,&redinitial,&redfinal,map);
     for (c=0;c<256;c++) scan->B[c] = redB[c][0];
     free (redB[0]);
     if (wlen > 0) /* reverse */
        {       /* reverse arrows */
          regularReverseArrows (redtrans,scan->dm,redinitial,redfinal,
                                &revtrans,&revinitial,&revfinal);
          free (redtrans[0]); free (redtrans);
          free (redinitial); free (redfinal);
                /* make all states initial */
          for (i=0;i<scan->dm;i++) SET(revinitial,i);
                /* make deterministic */
          regularMakeDet1 (scan->width,revtrans,scan->dm,&scan->dtrans);
          scan->dinitial = revinitial[0]; scan->dfinal = revfinal[0];
          free (revtrans[0]); free (revtrans);
          free (revinitial); free (revfinal);
		/* load table for skip loop */
          this = (scan->width == W) ? ONES : (ONE << scan->width)-1;
	  for (c=0;c<256;c++)
	     { tmp = scan->dinitial & scan->B[c];
	       scan->dtransin[c] = ZEROS;
	       for (i=0;i<scan->slices;i++)
	          { scan->dtransin[c] |= scan->dtrans[i][tmp&this];
		    tmp >>= scan->width;
	          }
	      }
        }
     else   /* forward */
        {       /* add self loop on initial states (in fact just 0) */
          SET(redtrans[0],0);
          for (c=0;c<256;c++) scan->B[c] |= ONE;
                /* make deterministic */
          regularMakeDet1 (scan->width,redtrans,scan->dm,&scan->dtrans);
          scan->dinitial = redinitial[0]; scan->dfinal = redfinal[0];
          free (redtrans[0]); free (redtrans);
          free (redinitial); free (redfinal);
        }
     return scan;
   }

regularData *regularPreproc (byte *pat, Tree *tree, Tree **pos)

   { regularData *P = malloc (sizeof(regularData));
     int slices,i,j,c,wlen,*map;
     Mask *trans,*rtrans; /* nondet trans */
     Mask active,initial,final;
     Mask S[256],A;

	/* allocate and load the masks */

     P->m = 0;
     regularLength(tree,&P->m);
     P->initial = createMask(P->m);
     P->current = createMask(P->m);
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
     regularFindBest (tree,trans,rtrans,P->m,P->B,P->initial,P->final,&wlen,
		      &active,&initial,&final,0);
     P->type = (wlen == 0) ? FWD:BWD;

	/* create verification automata */
     slices = (P->m+OptDetWidth-1)/OptDetWidth;
     P->width = (P->m+slices-1)/slices;
     regularMakeDet (P->width,trans,P->m,&P->dftransV);
     regularMakeDet (P->width,rtrans,P->m,&P->dbtransV);

	/* create searching subexpression according to the type of search */
     switch (detClass(tree,active,P->m))
	{ case SIMPLE:
	     for (i=0;i<P->m;i++) if (ISSET(active,i)) break;
	     for (j=i;j<P->m;j++) if (!ISSET(active,j)) break;
	     P->scanData = simpleLoadFast(wlen?true:false,P->B,i,j);
	     if (wlen) P->match = ONE<<(i+1); else P->match = ONE<<j;
	     P->scanText = (scanProc)simpleScan;
	     P->scanFree = (freeScanProc)simpleFreeScan;
	     break;
	  case EXTENDED:
	     for (i=0;i<P->m;i++) if (ISSET(active,i)) break;
	     for (j=i;j<P->m;j++) if (!ISSET(active,j)) break;
	     P->scanData = extendedLoadFast(wlen,P->B,S,A,i,j);
	     if (wlen) P->match = ONE<<(i+1); else P->match = ONE<<j;
	     P->scanText = (scanProc)extendedScan;
	     P->scanFree = (freeScanProc)extendedFreeScan;
	     break;
	  case REGULAR:
	     P->scanData = regularLoadFast(P->m,trans,wlen,P->B,
					   active,initial,final,&map);
	     P->revMap = malloc(P->m*sizeof(int));
	     for (i=0;i<P->m;i++)
	         if (ISSET(active,i)) P->revMap[map[i]] = i;
	     free (map);
	     P->dm = ((regularScanData*)P->scanData)->dm;
	     P->scanText = (scanProc)regularScan;
	     P->scanFree = (freeScanProc)regularFreeScan;
	     break;
        }
     P->match = ZEROS;

     free (initial); free (final); free (active);
     free (trans[0]); free (trans);
     free (rtrans[0]); free (rtrans);
     P->V1 = createMask (P->m);
     return P;
   }

void regularFree (regularData *P)

	/* Frees P */

   { P->scanFree (P->scanData);
     free (P->revMap);
     free (P->B[0]);
     free (P->dftransV[0][0]); free (P->dftransV[0]); free (P->dftransV);
     free (P->initial); free (P->final);
     free (P->dbtransV[0][0]); free (P->dbtransV[0]); free (P->dbtransV);
     free (P->V1);
     free (P);
   }

static byte *fwdCheck (byte *ptr, byte *top, Mask B[256],
		       Mask **dtrans, Mask current, Mask final,
		       int m, int width, Mask tmp)

		/* forward check.
		   returns NULL if not found, last pos read if found */
		/* assumes that ptr needs to be shifted before reading */

   { int i,f;
     int slices = (m+width-1)/width;
     COPY(tmp,current,m);
     while (true)
	{ if (!ISZERO(AND(tmp,final,m),m)) /* match, check context */
	     { if (recCheckRightContext(ptr+1,top+1)) return ptr;
	     }
	  if (ptr == top) return NULL;  /* end of buffer */
	   /* now compute new state in tmp */
	  ZERO (tmp,m);
	  f = 0;
	  for (i=0;i<slices;i++)
	      { OR(tmp,dtrans[i][SLICE(current,f,width)],m);
	        f += width;
	      }
	  ptr++; AND(tmp,B[*ptr],m);
	  if (ISZERO(tmp,m)) return NULL;  /* automaton died */
	  COPY(current,tmp,m);
	}
	   
   }
static byte *bwdCheck (byte *ptr, byte *top, Mask *B,
		       Mask **dtrans, Mask current, Mask final,
		       int m, int width, Mask tmp)

		/* directional check, forward or backward.
		   returns NULL if not found, last pos read if found */
		/* assumes that ptr needs to be shifted before reading */

   { int i,f;
     int slices = (m+width-1)/width;
     COPY(tmp,current,m);
     while (true)
	{ if (!ISZERO(AND(tmp,final,m),m)) /* match, check context */
	     { if (recCheckLeftContext(ptr,top)) return ptr;
	     }
	  if (ptr == top) return NULL;  /* end of buffer */
	   /* now compute new state in tmp */
	  ptr--; AND(current,B[*ptr],m);
	  ZERO (tmp,m);
	  f = 0;
	  for (i=0;i<slices;i++)
	      { OR(tmp,dtrans[i][SLICE(current,f,width)],m);
	        f += width;
	      }
	  if (ISZERO(tmp,m)) return NULL;  /* automaton died */
	  COPY(current,tmp,m);
	}
	   
   }

static bool checkMatch (regularData *P, byte *pos, byte **beg, byte **end)

	/* Checks that an initially promising match starting at pos
	   (if P->wlen) or ending at pos-1 (if !P->wlen)
	   corresponds or not to a complete match. The buffer limits
	   are *beg and *end, and in case of successful match we return
	   there the boundaries of the match found. */

   { byte *rbeg,*rend,*obeg,*oend;
     int i;

        /* start by knowing my surrounding record */

     if (P->type == FWD)  /* pos-1 is last position read */
        { recGetRecord (pos-1,*beg,*end,&rbeg,&rend,&obeg,&oend);
          if ((rbeg > pos-1) || (rend <= pos-1)) return false; /* overlaps */
        }
     else
        { recGetRecord (pos,*beg,*end,&rbeg,&rend,&obeg,&oend);
          if ((rbeg > pos) || (rend <= pos)) return false; /* overlaps */
        }

	/* now check which states have matched and test one by one */
     for (i=0;i<P->dm;i++)
	 if (P->match & (ONE<<i))
	    if (P->type == FWD)
	       {   /* first part: backward check */
		 ZERO(P->current,P->m);
		 SET(P->current,P->revMap[i]);
		 if (bwdCheck (pos,rbeg,P->B,P->dbtransV,P->current,
			       P->initial,P->m,P->width,P->V1) == NULL)
		    continue;
	           /* second part: forward check */
		 ZERO(P->current,P->m);
		 SET(P->current,P->revMap[i]);
                 if (fwdCheck (pos-1,rend-1,P->B,P->dftransV,P->current,
			       P->final,P->m,P->width,P->V1) == NULL)
		    continue;
	           /* match found */
                 *beg = obeg;
                 *end = oend;
                 return true;
	       }
            else
	       {   /* first part: forward check */
		 ZERO(P->current,P->m);
		 SET(P->current,P->revMap[i]);
		 if (fwdCheck (pos,rend-1,P->B,P->dftransV,P->current,
			       P->final,P->m,P->width,P->V1) == NULL)
		    continue;
	           /* second part: backward check */
		 ZERO(P->current,P->m);
		 SET(P->current,P->revMap[i]);
                 if (bwdCheck (pos+1,rbeg,P->B,P->dbtransV,P->current,
			       P->initial,P->m,P->width,P->V1) == NULL)
		    continue;
	           /* match found */
                 *beg = obeg;
                 *end = oend;
                 return true;
	       }
     return false;
   }

static bool fwdScan1 (byte **beg, byte **end, checkProc checkMatch, 
		      regularData *P, mask **dtrans, register mask *B,
		      mask initial, register mask final)
 
   { register mask current;
     register byte *pos,*top;	/* current and final text pointers */
     register mask *dtrans0 = dtrans[0];

     pos = *beg; top = *end;
     while (pos < top)
        { current = initial;
          while (pos < top)
	     { if (current & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = current & final;
	            if (checkMatch (P,pos,beg,end)) return true;
                  }
		     /* take new character and compute transition */
	       current = dtrans0[current] & B[*pos++];
	     }
	}
     return false;
   }

static bool fwdScan2 (byte **beg, byte **end, checkProc checkMatch, 
		      regularData *P, mask **dtrans, register mask *B,
		      mask initial, register mask final, register int width)
 
   { register mask current;
     register byte *pos,*top;	/* current and final text pointers */
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask *dtrans0 = dtrans[0];
     register mask *dtrans1 = dtrans[1];

     pos = *beg; top = *end;
     while (pos < top)
        { current = initial;
          while (pos < top)
	     { if (current & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = current & final;
	            if (checkMatch (P,pos,beg,end)) return true;
                  }
		     /* compute transition */
	       current = (dtrans0[current&this]|dtrans1[current>>width]) & B[*pos++];
	     }
	}
     return false;
   }

static bool fwdScan (byte **beg, byte **end, checkProc checkMatch, 
		     regularData *P, register mask **dtrans, register mask *B,
		     mask initial, register mask final,
		     register int slices, register int width)
 
   { register mask current,tmp;
     register byte *pos,*top;	/* current and final text pointers */
     register int c,i;
     register mask this = (width == W) ? ONES : (ONE << width)-1;

     pos = *beg; top = *end;
     while (pos < top)
        { current = initial;
          while (pos < top)
	     { if (current & final)
	          {  /* selected part of pattern has matched, check all */
		    P->match = current & final;
	            if (checkMatch (P,pos,beg,end)) return true;
                  }
	       if ((c = *pos++) == OptRecChar) break; // faster
		     /* compute transition */
	       tmp = current; current = ZEROS;
	       for (i=0;i<slices;i++)
	           { current |= dtrans[i][tmp&this];
		     tmp >>= width;
	           }
	       current &= B[c];
	     }
	}
     return false;
   }

static bool bwdScan1 (byte **beg, byte **end, checkProc checkMatch, 
		      regularData *P, mask **dtrans, register mask *B,
		      register mask final, register mask *dtransin,
		      register int wlen)
 
   { register mask icurrent,current;
     register byte *pos,*top;	/* current and final text pointers */
     register mask *dtrans0 = dtrans[0];
     register int j;

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
			 if (checkMatch (P,pos+1,beg,end)) return true;
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

static bool bwdScan2 (byte **beg, byte **end, checkProc checkMatch, 
		      regularData *P, mask **dtrans, register mask *B,
		      register mask final, register mask *dtransin,
		      register int width, register int wlen)
 
   { register mask current,icurrent;
     register byte *pos,*top;	/* current and final text pointers */
     register mask this = (width == W) ? ONES : (ONE << width)-1;
     register mask *dtrans0 = dtrans[0];
     register mask *dtrans1 = dtrans[1];
     register int j;

     pos = *beg; top = *end;
     pos--;
     top -= wlen;
     while (pos < top)
        { if (!(current = dtransin[pos[wlen]]))  /* skip loop */
	     { pos += wlen; continue; }
	  j = wlen-1;
	  while (true)
	     {   /* compute transition */
	       icurrent = current & B[pos[j]];
	       current = dtrans0[icurrent&this]|dtrans1[icurrent>>width];
	       if (!current) break;
	       if (--j == 0) 
	          {  /* a factor of the pattern has matched, 
			if it is a prefix, check all */
	            if (current & final)
		       { P->match = icurrent;
			 if (checkMatch (P,pos+1,beg,end)) return true;
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

static bool bwdScan (byte **beg, byte **end, checkProc checkMatch, 
		     regularData *P, register mask **dtrans, register mask *B,
		     register mask final, register mask *dtransin,
		     register int slices, register int width, register int wlen)
 
   { register mask current,tmp,icurrent;
     register byte *pos,*top;	/* current and final text pointers */
     register int ch,i,j;
     register mask this = (width == W) ? ONES : (ONE << width)-1;

     pos = *beg; top = *end;
     pos--;
     top -= wlen;
     while (pos < top)
        { if (!(current = dtransin[pos[wlen]]))  /* skip loop */
	     { pos += wlen; continue; }
	  j = wlen-1;
	  while (true)
	     {   /* compute transition */
	       icurrent = tmp = current & B[pos[j]]; current = ZEROS;
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
			 if (checkMatch (P,pos+1,beg,end)) return true;
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

bool regularScan (byte **beg, byte **end, checkProc checkMatch, void *P,
		  regularScanData *scan)
 
        /* Scans from *beg to *end-1 for P. Returns if it could find
           it. In case of returning true, *beg and *end are set to limit
           the first occurrence of P. It requires
             - a procedure
                bool checkMatch (P,byte *pos, byte **beg, byte **end)
               which whether there is a match around pos and returns in
               *beg and *end the limits. "Around" means that pos is the initial
	       position of the search pattern given to regularScan.
             -  dtrans, initial, final, m, width: data for the automaton
	     - wlen (zero => fwd scan)
         */                                                                    

   { if (scan->wlen > 0)  /* bwd scan */
	{ if (scan->slices == 1)
	     return bwdScan1 (beg,end,checkMatch,P,scan->dtrans,scan->B,
			      scan->dfinal,scan->dtransin,scan->wlen);
	  if (scan->slices == 2)
	     return bwdScan2 (beg,end,checkMatch,P,scan->dtrans,scan->B,
			    scan->dfinal,scan->dtransin,scan->width,scan->wlen);
	  return bwdScan (beg,end,checkMatch,P,scan->dtrans,scan->B,
			  scan->dfinal,scan->dtransin,scan->slices,scan->width,
			  scan->wlen);
	}
     else  /* fwd scan */
	{ if (scan->slices == 1)
	     return fwdScan1 (beg,end,checkMatch,P,scan->dtrans,scan->B,
			      scan->dinitial,scan->dfinal);
	  if (scan->slices == 2)
	     return fwdScan2 (beg,end,checkMatch,P,scan->dtrans,scan->B,
			      scan->dinitial,scan->dfinal,scan->width);
	  return fwdScan (beg,end,checkMatch,P,scan->dtrans,scan->B,
			  scan->dinitial,scan->dfinal,scan->slices,scan->width);
	}
   }

bool regularSearch (byte **beg, byte **end, regularData *P)

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the first occurrence of P */

   { return P->scanText (beg,end,(checkProc)checkMatch,P,P->scanData);
   }
 
