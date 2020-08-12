UNAME := $(shell uname)
LBITS := $(shell getconf LONG_BIT)
ifeq ($(UNAME), Darwin)
        CC=`which clang`
        LIB=../libmac
else ifeq ($(UNAME),Linux)
        CC=`which gcc`
        ifeq ($(LBITS),64)
                #LIB=../lib64
        else
                LIB=../lib32
        endif
else
	UNAME := $(shell uname -o) 
	ifeq ($(UNAME), Cygwin)
		CC=`which gcc`
		ifeq ($(LBITS),64)
			LIB=../libcygwin64
		else
			LIB=../libcygwin32
		endif
	endif
endif

all:	tephra2-2012
		mv tephra2-2012 ../tephra2-2012

tephra2-2012:	new_tephra.o ../common_src/parameters.h tephra2_calc.o
		$(CC) -Wall -O3 -idirafter ../include -o tephra2-2012 new_tephra.o tephra2_calc.o -lm $(LIB)/libgc.a -ldl

new_tephra.o:	new_tephra.c ../common_src/common_structures.h makefile
		$(CC) -Wall -O3 -idirafter ../include -c new_tephra.c

tephra2_calc.o:		../common_src/tephra2_calc.c ../common_src/common_structures.h makefile
		$(CC) -Wall -O3 -idirafter ../include -c ../common_src/tephra2_calc.c

clean:
		rm -fv *.o
