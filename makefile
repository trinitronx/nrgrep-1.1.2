#Nrgrep -- a fast and flexible pattern matching tool.
#Copyright (C) 2000 Gonzalo Navarro
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of
#Chile. Blanco Encalada 2120, Santiago, Chile. gnavarro@dcc.uchile.cl

FLAGS = -O9

nrgrep: shell.o record.o search.o eregular.o regular.o eextended.o extended.o esimple.o simple.o parser.o bitmasks.o options.o memio.o buffer.o makefile
	gcc $(FLAGS) -o nrgrep shell.o record.o search.o eregular.o regular.o eextended.o extended.o esimple.o simple.o parser.o bitmasks.o options.o memio.o buffer.o
	strip nrgrep

buffer.o: buffer.c buffer.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c buffer.c

memio.o: memio.c memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c memio.c

options.o: options.c options.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c options.c

bitmasks.o: bitmasks.c bitmasks.h basics.h memio.h makefile
	gcc $(FLAGS) -c bitmasks.c

parser.o: parser.c record.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c parser.c

simple.o: simple.c simple.h record.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c simple.c

esimple.o: esimple.c esimple.h simple.h record.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c esimple.c

extended.o: extended.c extended.h simple.h record.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c extended.c

eextended.o: eextended.c eextended.h extended.h esimple.h record.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c eextended.c

regular.o: regular.c regular.h extended.h simple.h record.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c regular.c

eregular.o: eregular.c eregular.h regular.h extended.h simple.h record.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c eregular.c

search.o: search.c search.h regular.h eextended.h extended.h esimple.h simple.h record.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c search.c

record.o: record.c record.h search.h simple.h buffer.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c record.c

shell.o: shell.c record.h search.h simple.h buffer.h parser.h options.h memio.h except.h bitmasks.h basics.h makefile
	gcc $(FLAGS) -c shell.c

