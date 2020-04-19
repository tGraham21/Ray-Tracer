CFLAGS  = -g -Wall
CC = g++ $(CFLAGS)

a1: assignment1.o vector.o
	$(CC) -o a1 assignment1.o

assignment1.o: assignment1.cc
	$(CC) -c assignment1.cc

vector.o: vector.cc a1_struct.h
	$(CC) -c vector.cc

clean:
	rm -f *.o
