CC :=`which mpicc`
#CC=/usr/bin/mpicc
UNAME := $(shell uname)
LBITS := $(shell getconf LONG_BIT)
ifeq ($(UNAME), Darwin)
        LIB=../libmac
else ifeq ($(UNAME),Linux)
        ifeq ($(LBITS),64)
                LIB=../lib64
        else
                LIB=../lib32
        endif
else
        UNAME := $(shell uname -o)
        ifeq ($(UNAME), Cygwin)
                ifeq ($(LBITS),64)
                        LIB=../libcygwin64
                else
                        LIB=../libcygwin32
                endif
        endif
endif

all:  	tephra2012_inversion
	mv tephra2012_inversion ..

tephra2012_inversion:		master.o slave.o ameoba.o tephra_inversion.o minimizing_func.o tephra2_calc.o fit_tests.o
		$(CC)  -O3 -Wall -o tephra2012_inversion \
		master.o\
		slave.o\
		ameoba.o\
		fit_tests.o -lm\
		minimizing_func.o -lm\
		tephra_inversion.o\
		tephra2_calc.o $(LIB)/libgc.a -ldl

master.o:		master.c ../common_src/parameters.h makefile
			$(CC) -Wall  -O3 -idirafter ../include -c master.c

slave.o:		slave.c ../common_src/parameters.h makefile
			$(CC) -Wall -O3 -idirafter ../include -c slave.c

ameoba.o:		ameoba.c ../common_src/parameters.h makefile
			$(CC) -Wall -idirafter ../include -c ameoba.c

fit_tests.o:		fit_tests.c ../common_src/common_structures.h makefile
			$(CC) -Wall -idirafter ../include -c fit_tests.c
			
minimizing_func.o:	minimizing_func.c ../common_src/common_structures.h makefile 
			$(CC) -Wall -idirafter ../include -c minimizing_func.c

tephra2_calc.o:	../common_src/tephra2_calc.c ../common_src/common_structures.h makefile
			$(CC) -Wall -O3 -idirafter ../include -c ../common_src/tephra2_calc.c

tephra_inversion.o:	tephra_inversion.c ../common_src/common_structures.h ../common_src/parameters.h ../common_src/prototypes.h makefile
			$(CC) -Wall -O3 -idirafter ../include -c tephra_inversion.c

clean:
	rm *.o
