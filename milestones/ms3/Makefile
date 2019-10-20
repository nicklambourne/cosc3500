CC = g++
MPICC = mpic++
CFLAGS = -g -Wall -Werror -std=c++11

all: clean serial hybrid

serial:
	${CC} ${CFLAGS} lcs-serial.cc -o ./bin/lcs-serial

hybrid:
	${MPICC} ${CFLAGS} lcs-hybrid.cc -lm -fopenmp -o ./bin/lcs-hybrid

.PHONY: clean
clean:
	rm -rf *.o bin
	mkdir bin
