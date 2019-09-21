CC = g++
CFLAGS = -g -Wall -Werror -std=c++11

all: clean serial omp

serial:
	${CC} ${CFLAGS} lcs-serial.cc -o ./bin/lcs-serial

omp:
	${CC} ${CFLAGS} lcs-omp.cc -fopenmp -o ./bin/lcs-omp 


.PHONY: clean
clean:
	rm -rf *.o bin
	mkdir bin
