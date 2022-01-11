
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

	/* Main shell */

#include "basics.h"
#include "except.h"
#include "buffer.h"
#include "parser.h"
#include "record.h"
#include "simple.h"
#include "search.h"
#ifdef LINUX
#include <getopt.h>
#endif

void usage (char **argv)

   { fprintf(stderr,"Usage: %s [-options] <pattern> <list of files>\n"
	            "       %s -H for help\n",argv[0],argv[0]);
     exit(1);
   }

void help (char **argv, FILE *f)

   { fprintf (f,"\n");
     fprintf (f,"Nrgrep -- a fast and flexible pattern matching tool\n",argv[0]);
     fprintf (f,"Nrgrep version 1.1, Copyright (C) 2000,2001 Gonzalo Navarro\n");
     fprintf (f,"Nrgrep comes with ABSOLUTELY NO WARRANTY and is free software,\n");
     fprintf (f,"you are welcome to redistribute it under certain conditions.\n");
     fprintf (f,"See the GNU General Public License for details.\n");
     fprintf (f,"Please send bug reports to bug-nrgrep@gnu.org.\n");
     fprintf (f,"\nUsage: %s [-iclGhnvdbmskL] <pattern> <list of files>\n",argv[0]);
     fprintf (f,"         searches <pattern> in the list of files\n");
     fprintf (f,"         printing the regions surrounding the occurrences\n");
     fprintf (f,"\nOptions:\n");
     fprintf (f,"     -i: the search is case insensitive\n");
     fprintf (f,"     -w: only matches whole words\n");
     fprintf (f,"     -x: only matches whole records\n");
     fprintf (f,"     -c: just counts the matches, does not print them\n");
     fprintf (f,"     -l: output filenames only, not their contents\n");
     fprintf (f,"     -G: output whole files\n");
     fprintf (f,"     -h: do not output file names\n");
     fprintf (f,"     -n: output records preceded by record number\n");
     fprintf (f,"     -v: report nonmatching records\n");
     fprintf (f,"     -d <delim>: sets the record delimiter to <delim>\n");
     fprintf (f,"          Default is \\n. Use a # at the end to make the\n");
     fprintf (f,"          delimiter part of previous record (default is next)\n");
     fprintf (f,"     -b <bufsize>: sets the buffer size to <bufsize> in Kb\n");
     fprintf (f,"           Default is %i\n",OptBufSize);
     fprintf (f,"     -m <bits>: sets the maximum table sizes to 2^<bits> words\n");
     fprintf (f,"           Default is %i\n",OptDetWidth);
     fprintf (f,"     -s <sep>: sets the output record separator to <sep>\n");
     fprintf (f,"     -k <err>[idst]: allow up to <err> errors in the matches\n");
     fprintf (f,"         [idst] means permitting ins, del, subs, transp operations\n");
     fprintf (f,"         (default is all)\n");
     fprintf (f,"     -L: take pattern literally (no metachars)\n");
     fprintf (f,"\nPattern syntax:\n");
     fprintf (f,"  In the simplest form, the pattern is a sequence of letters.\n");
     fprintf (f,"  In particular, you can use '\\t' as a tab, '\\n' as a newline, '\\xdd' to\n");
     fprintf (f,"  denote the character with hex ASCII code dd, and '\\c' to make nrgrep\n");
     fprintf (f,"  interpret 'c' literally (e.g. '\\\\').\n");
     fprintf (f,"  In addition, a set of characters can be specified by enclosing it in square\n");
     fprintf (f,"  brackets (e.g. '[abc]' matches 'a', 'b' or 'c'). If '^' is put as the first\n");
     fprintf (f,"  character of the class, the class is complemented. A range of characters\n");
     fprintf (f,"  (e.g. 'A-Z') can be put inside a class. Other possible classes that are\n");
     fprintf (f,"  expressed without brackets are '#' (any separator) and '.' (any character).\n");
     fprintf (f,"  These simple patterns can be extended by appending operators to single\n");
     fprintf (f,"  letters or classes: '?' means that the preceding character is optional,\n");
     fprintf (f,"  '*' that it can appear zero or more times, and '+' that it can appear one\n");
     fprintf (f,"  or more times.\n");
     fprintf (f,"  Finally, the most complex patterns that can be expressed are the regular\n");
     fprintf (f,"  expressions, which also permit the union operator '|' (e.g. 'abc|de' matches\n");
     fprintf (f,"  the strings 'abc' and 'de') and the parenthesis '(' ')' to enclose\n");
     fprintf (f,"  subexpressions, so that '?', '*' and '+' can be applied to complete\n");
     fprintf (f,"  expressions and not only letters, e.g. 'ab(cd|e)*fg?h'.\n\n");
   }

main (int argc, char **argv)

   { char **opts;
     byte *patt;
     char **files;
     int f;
     Buffer B;
     searchData *sData;
     int matches;
     int opt;
     int i,j;

	/* get the options */

     opterr = 0;
     while ((opt = getopt(argc,argv,"HiwxclGhnvd:s:b:m:k:L")) != -1)
          switch (opt)
		 { case '?':     /* help and error */
		      help(argv,stderr);
		      error1 ("Unrecognized option %c",optopt);
		      break;
		   case 'H':    /* help and exit */
		      help(argv,stdout);
		      exit(0);
		   case 'i':    /* case insensitive search */
		      OptCaseInsensitive = true;
		      break;
		   case 'w':    /* whole word matching */
		      OptWholeWord = true;
		      break;
		   case 'x':    /* whole record matching */
		      OptWholeRecord = true;
		      break;
		   case 'c':    /* just count the matches */
		      OptRecPrint = false;
		      break;
		   case 'l':    /* output filenames only */
		      OptRecFiles = true;
		      break;
		   case 'G':    /* output whole files only */
		      OptRecPrintFiles = true;
		      break;
		   case 'n':    /* output records preceded by record number */
		      OptRecNumber = true;
		      break;
		   case 'h':    /* do not output file names */
		      OptRecFileNames = false;
		      break;
		   case 'v':    /* report nonmatching records */
		      OptRecPositive = false;
		      break;
		   case 'd':	/* record delimiter */
		      if (optarg == NULL) /* (argc < 2) */
			 error0 ("A delimiter expected after -d");
		      OptRecPatt = optarg; /* opts[1]; */
		      OptRecPos = 0;
		      if (OptRecPatt[strlen(OptRecPatt)-1] == '#')
			 { OptRecPos = strlen(OptRecPatt)-1;
			   OptRecPatt[OptRecPos] = 0;
			 }
		      if (strlen(OptRecPatt) == 1)
			 OptRecChar = OptRecPatt[0];
		      else OptRecChar = -1;
		      break;
		   case 's':	/* output record separator */
		      if (optarg == NULL) /* (argc < 2) */
			 error0 ("A separator expected after -s");
		      OptRecSep = malloc (strlen(optarg));
		      i = 0; j = 0;
		      while (optarg[i]) /* opts[1] */
			 OptRecSep[j++] = getAchar (optarg,&i);
		      OptRecSep[j] = 0;
		      break;
		   case 'b':	/* buffer size */
		      if (optarg == NULL) /* (argc < 2) */
			 error0 ("A buffer size expected after -b");
		      OptBufSize = atoi(optarg); /* atoi(opts[1]); */
		      break;
		   case 'm':	/* tables width */
		      if (optarg == NULL) /* (argc < 2) */
			 error0 ("A number of bits expected after -m");
		      i = atoi(optarg); /* atoi(opts[1]); */
		      if ((i<=0) || (i>W) || (W % i))
			 { warn2("The number of bits must be between 1 and %i "
				"and divide %i, after -m",W,W);
			 }
		      else OptDetWidth = i;
		      break;
		   case 'k':	/* errors */
		      if (optarg == NULL) /* (argc < 2) */
			 error0 ("A number of errors expected after -k");
		      if (optarg[0] && !isdigit(optarg[strlen(optarg)-1]))
			 { OptIns = OptDel = OptSubs = OptTransp = false;
		           do { switch (optarg[strlen(optarg)-1])
				   { case 'i': OptIns = true; break;
				     case 'd': OptDel = true; break;
				     case 's': OptSubs = true; break;
				     case 't': OptTransp = true; break;
				     default: error0 ("<num>[idst] expected after -k");
			           }
				optarg[strlen(optarg)-1] = 0;
			      }
		           while (optarg[0] && !isdigit(optarg[strlen(optarg)-1]));
			 }
		      OptErrors = atoi(optarg); /* atoi(opts[1]); */
		      break;
		   case 'L':    /* take pattern literally */
		      OptLiteral = true;
		      break;
		}
     opts = argv+optind; argc -= optind;

		/* some consistency checks */
     if ((argc <= 2) && !OptRecFiles) OptRecFileNames = false;
     if (!OptRecPrint && OptRecPrintFiles)
	{ warn0 ("Options -c and -G are incompatible, overriding -G");
	  OptRecPrintFiles = false;
	}
     if (OptRecNumber && OptRecFiles)
	{ warn0 ("Options -n and -l are incompatible, overriding -n");
	  OptRecNumber = false;
	}
     if (OptRecFiles && OptRecPrintFiles)
	{ warn0 ("Options -l and -G are incompatible, overriding -l");
	  OptRecFiles = false;
	}
	
	/* get the pattern */

     if (argc < 1) usage (argv);
     patt = opts[0];
     if (patt[0] == '^') 
        { OptStartLine = true; patt++; }
     if (patt[strlen(patt)-1] == '$') 
        { OptEndLine = true; patt[strlen(patt)-1] = 0; }
     sData = searchPreproc (patt);
     if (sData == NULL) error1("Syntax error in pattern %s",patt);
     recPreproc ();

	/* search the files */

     files = opts+1; argc--;
     B = bufCreate();
     matches = 0;
     if (argc == 0)  /* stdin */
	{ if (OptRecPrintFiles)
	     { warn0 ("Option -G does not work on standard input, overriding it");
               OptRecPrintFiles = false;
	     }
	  f = fileno (stdin);
	  bufSetFile (B,f);
	  if (OptRecNumber || !OptRecPositive)
	    matches += recSearchRecFile (NULL,B,sData);
	  else matches += recSearchFile (NULL,B,sData);
	}
     else
        while (argc)
	   { int f = open (*files,O_RDONLY);
	     int newm;
	     if (f == -1) error1 ("Cannot read file %s",*files);
	     bufSetFile (B,f);
	     if (OptRecNumber || !OptRecPositive)
	          newm = recSearchRecFile (*files,B,sData);
	     else newm = recSearchFile (*files,B,sData);
	     close (f);
	     matches += newm;
	     if (OptRecPrintFiles && (OptRecPositive == (newm != 0)))
		{ char str[1024];
		  sprintf (str,"cat %s",*files);
		  system (str);
		  printf (OptRecSep);
		}
	     files++; argc--;
	   }

	/* final report */

     if ((OptRecFiles && !OptRecFileNames) ||
	 (!OptRecFiles && !OptRecPrint)) /* report # of occs */
	printf ("Total: %i matching %s\n", matches,
		OptRecFiles ? "files" : "records");

	/* clean up */

     bufDestroy (B);
     searchFree (sData);
     recFree ();
     exit(0);
   } 
