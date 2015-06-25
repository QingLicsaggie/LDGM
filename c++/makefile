CC = g++ -std=c++0x
all: ldgm
ldgm: ldgm.o main.o
	$(CC) -o ldgm ldgm.o main.o
ldgm.o: ldgm.hpp ldgm.cpp
	$(CC) -c ldgm.cpp
main.o: main.cpp ldgm.cpp ldgm.hpp
	$(CC) -c main.cpp
clean:
	/bin/rm -f *.o ldgm