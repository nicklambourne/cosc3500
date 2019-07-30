CXX = gcc                                       
CFLAGS = -g -Wall -Werror -pedantic

all: clean main postclean

main: main.o
	$(CXX) $(CFLAGS) main.o -o main

main.o: main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

.PHONY: clean
clean:
	rm -f *.o main

.PHONY: postclean
postclean:
	rm -f *.o