CC=/usr/bin/mpicc
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
   # do 64 bit stuff here, like set some CFLAGS
	LIB=../lib64
else
   # do 32 bit stuff here
	LIB=../lib32
endif

all:  	tephra2012_inversion
	mv tephra2012_inversion ..

tephra2012_inversion:		master.o slave.o ameoba.o tephra_inversion.o minimizing_func.o tephra2_calc.o
		$(CC)  -O3 -Wall -pthread -o tephra2012_inversion \
		master.o\
		slave.o\
		ameoba.o\
		minimizing_func.o -lm\
		tephra_inversion.o\
		tephra2_calc.o $(LIB)/gc.a -ldl

master.o:		master.c ../common_src/parameters.h makefile
			$(CC) -Wall  -O3 -idirafter ../include -c master.c

slave.o:		slave.c ../common_src/parameters.h makefile
			$(CC) -Wall -O3 -idirafter ../include -c slave.c

ameoba.o:		ameoba.c ../common_src/parameters.h makefile
			$(CC) -Wall -idirafter ../include -c ameoba.c

minimizing_func.o:	minimizing_func.c ../common_src/common_structures.h makefile 
			$(CC) -Wall -idirafter ../include -c minimizing_func.c

tephra2_calc.o:	../common_src/tephra2_calc.c ../common_src/common_structures.h makefile
			$(CC) -Wall -O3 -idirafter ../include -c ../common_src/tephra2_calc.c

tephra_inversion.o:	tephra_inversion.c ../common_src/common_structures.h ../common_src/parameters.h ../common_src/prototypes.h makefile
			$(CC) -Wall -O3 -idirafter ../include -c tephra_inversion.c

clean:
	rm *.o tephra2012_inversion
