all: clean serial omp

serial:
	g++ -g -Wall -Werror -std=c++11 lcs-serial.cc -o ./bin/lcs-serial

omp:
	g++ −g −Wall -Werror -std=c++11 lcs-omp.cc −fopenmp −o ./bin/lcs-omp

.PHONY: clean
clean:
	rm -rf *.o bin
	mkdir bin
