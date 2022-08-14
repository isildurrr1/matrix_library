default: all

all: clean build

build: matrix.a

matrix.a: matrix.o
	ar rcs matrix.a matrix.o
	ranlib matrix.a

matrix.o:  matrix.c matrix.h
	gcc -c matrix.c

clean:
	-rm -rf *.o *.a

rebuild: clean all

valgrind:
	valgrind --leak-check=full -s -q ./test

