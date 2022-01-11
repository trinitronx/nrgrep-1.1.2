
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

	/* Search a pattern in a text, general interface */

#include "search.h"
#include "record.h"
#include "simple.h"
#include "esimple.h"
#include "extended.h"
#include "eextended.h"
#include "regular.h"
#include "eregular.h"
#include "parser.h"

searchData *searchPreproc (byte *pat)

	/* Preprocesses pat and creates a searchData structure
	   for searching */

   { searchData *P;
     Tree *tree;
     Tree **pos;
     int m,class;

     m = strlen(pat);
     pos = malloc (m * sizeof(Tree*));
     tree = parse (pat, m, pos);
     if ((tree == NULL) || tree->eps)
	{ freeTree(tree); free (pos); return NULL; }
     class = detClass (tree,NULL,0);

     P = malloc (sizeof(searchData));

     switch (class)
        { case SIMPLE:
            if (!OptErrors)
               { P->preprocData = simplePreproc (pat,tree,pos);
	         P->searchType = SIMPLE;
	         P->search = simpleSearch;
	       }
            else
               { P->preprocData = esimplePreproc (pat,tree,pos,OptErrors);
	         P->searchType = ESIMPLE;
	         P->search = esimpleSearch;
	       }
	    break;
          case EXTENDED:
            if (!OptErrors)
               { P->preprocData = extendedPreproc (pat,tree,pos);
	         P->searchType = EXTENDED;
	         P->search = extendedSearch;
	       }
            else
               { P->preprocData = eextendedPreproc (pat,tree,pos,OptErrors);
	         P->searchType = EEXTENDED;
	         P->search = eextendedSearch;
	       }
	    break;
          case REGULAR:
            if (!OptErrors)
               { P->preprocData = regularPreproc (pat,tree,pos);
	         P->searchType = REGULAR;
	         P->search = regularSearch;
	       }
            else
               { P->preprocData = eregularPreproc (pat,tree,pos,OptErrors);
	         P->searchType = EREGULAR;
	         P->search = eregularSearch;
	       }
	    break;
        }
     free (pos); freeTree (tree);
     return P;
   }

void searchFree (searchData *P)

	/* Frees P */

   { switch (P->searchType)
	{ case SIMPLE:
		simpleFree (P->preprocData);
		break;
	  case ESIMPLE:
		esimpleFree (P->preprocData);
		break;
	  case EXTENDED:
		extendedFree (P->preprocData);
		break;
	  case EEXTENDED:
		eextendedFree (P->preprocData);
		break;
  	  case REGULAR:
		regularFree (P->preprocData);
		break;
  	  case EREGULAR:
		eregularFree (P->preprocData);
		break;
        }
     free (P);
   }
		
bool searchScan (byte **beg, byte **end, searchData *P)

	/* Searches from *beg to *end-1 for P. Returns if it could find
	   it. In case of returning true, *beg and *end are set to limit
	   the printable record containing the occurrence of P */

   { return P->search (beg,end,P->preprocData);
   }

