all: main_j.out

main_j.out: main_j.o fonctions_j.o
	gcc -Wall -g -o main_j.out main_j.o fonctions_j.o

main_j.o: main_j.c fonctions_j.h
	gcc -Wall -g -c -o main_j.o main_j.c

fonctions_j.o: fonctions_j.c fonctions_j.h
	gcc -Wall -g -c -o fonctions_j.o fonctions_j.c

clean:
	rm *.o *.out

mrproper: clean
	rm -rf main_j.out
