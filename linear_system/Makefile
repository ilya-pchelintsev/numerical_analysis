all: solver

solver: main.o matrread.o linsolve.o utils.o
	gcc -O3 main.o matrread.o linsolve.o utils.o -o solver

main.o: main.c
	gcc -O3 -c main.c

matrread.o: matrread.c
	gcc -O3 -c matrread.c

linsolve.o: linsolve.c
	gcc -O3 -c linsolve.c

utils.o: utils.c
	gcc -O3 -c utils.c

clean:
	rm -rf *.o solver

