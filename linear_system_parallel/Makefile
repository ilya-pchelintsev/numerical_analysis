all: solver

FLAGS=

solver: main.o matrread.o linsolve.o utils.o
	gcc $(FLAGS) main.o matrread.o linsolve.o utils.o -o solver

main.o: main.c
	gcc $(FLAGS) -c main.c

matrread.o: matrread.c
	gcc $(FLAGS) -c matrread.c

linsolve.o: linsolve.c
	gcc $(FLAGS) -c linsolve.c

utils.o: utils.c
	gcc $(FLAGS) -c utils.c

clean:
	rm -rf *.o solver

