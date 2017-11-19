all: main.out

main.out: main.o fonctions.o
	gcc -Wall -o main.out main.o fonctions.o

main.o: main.c fonctions.h
	gcc -Wall -c -o main.o main.c

fonctions.o: fonctions.c fonctions.h 
	gcc -Wall -c -o fonctions.o fonctions.c

clean:
	rm *.o *.out

mrproper: clean
	rm -rf main.out
