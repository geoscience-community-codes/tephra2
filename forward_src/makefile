CC=/usr/bin/gcc
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
   # do 64 bit stuff here, like set some CFLAGS
   LIB=../lib64
else
   # do 32 bit stuff here
   LIB=../lib32
endif

all:	tephra2-2012
		mv tephra2-2012 ../tephra2-2012

tephra2-2012:	new_tephra.o ../common_src/parameters.h tephra2_calc.o
		$(CC) -Wall -O3 -pthread -idirafter ../include -o tephra2-2012 new_tephra.o tephra2_calc.o -lm $(LIB)/gc.a -ldl

new_tephra.o:	new_tephra.c ../common_src/common_structures.h makefile
		$(CC) -Wall -O3 -idirafter ../include -c new_tephra.c

tephra2_calc.o:		../common_src/tephra2_calc.c ../common_src/common_structures.h makefile
		$(CC) -Wall -O3 -idirafter ../include -c ../common_src/tephra2_calc.c

clean:
		rm -fv *.o ../tephra2-2012
