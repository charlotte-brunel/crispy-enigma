all: main_j.out

main_j.out: main_j.o fonctions_j.o
	gcc -Wall -o main_j.out main_j.o fonctions_j.o

main_j.o: main_j.c fonctions_j.h
	gcc -Wall -c -o main_j.o main_j.c

fonctions_j.o: fonctions_j.c fonctions_j.h
	gcc -Wall -c -o fonctions_j.o fonctions_j.c

clean:
	rm *.o *.out

mrproper: clean
	rm -rf main_j.out
