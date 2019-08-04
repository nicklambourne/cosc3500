all: clean bin/suffix_tree

bin/suffix_tree:
	g++ -Wall -Werror -std=c++11 suffix_tree.cpp -o ./bin/suffix_tree

.PHONY: clean
clean:
	rm -rf *.o bin
	mkdir bin
