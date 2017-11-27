all: main.out

main.out: main.o function.o
	gcc -Wall -o main.out main.o function.o

main.o: main.c function.h
	gcc -Wall -c -o main.o main.c

function.o: function.c function.h
	gcc -Wall -c -o function.o function.c

clean:
	rm *.o *.out

mrproper: clean
rm -rf main.out
