all: clean bin/lcs

bin/lcs:
	g++ -g -Wall -Werror -std=c++11 main.cpp -o ./bin/lcs

.PHONY: clean
clean:
	rm -rf *.o bin
	mkdir bin
