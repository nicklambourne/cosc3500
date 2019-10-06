CC = g++
MPICC = mpic++
CFLAGS = -g -Wall -Werror -std=c++11

all: clean serial omp hybrid

serial:
	${CC} ${CFLAGS} lcs-serial.cc -o ./bin/lcs-serial

omp:
	${CC} ${CFLAGS} lcs-omp.cc -fopenmp -o ./bin/lcs-omp 

hybrid:
	${MPICC} ${CFLAGS} lcs-hybrid.cc -lm -O3 -fopenmp -o ./bin/lcs-hybrid

.PHONY: clean
clean:
	rm -rf *.o bin
	mkdir bin
