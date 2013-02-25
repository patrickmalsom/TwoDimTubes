CC=gcc
#order of arguments: Temp, sigma, sigma^2
CFLAGS=-Wall -Warray-bounds -fopenmp -march=native -O2 -g
LFLAGS=-fopenmp -lm -lgsl -lgslcblas

all: TwoDimTubes

TwoDimTubes.o: TwoDimTubes.c
	${CC} ${CFLAGS} ${EXTFLAGS} -c TwoDimTubes.c

TwoDimTubes: TwoDimTubes.o 
	${CC} ${LFLAGS} TwoDimTubes.o -o TwoDimTubes.out
 
clean: 
	rm -f TwoDimTubes.o TwoDimTubes.out

