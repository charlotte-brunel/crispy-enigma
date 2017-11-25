all: main_algo_test.out

main_algo_test.out: main_algo_test.o fonctions_algo_test.o
	gcc -Wall -o main_algo_test.out main.o fonctions_algo_test.o

main_algo_test.o: main_algo_test.c fonctions_algo_test.h
	gcc -Wall -c -o main_algo_test.o main_algo_test.c

fonctions_algo_test.o: fonctions_algo_test.c fonctions_algo_test.h
	gcc -Wall -c -o fonctions_algo_test.o fonctions_algo_test.c

clean:
	rm *.o *.out

mrproper: clean
	rm -rf main_algo_test.out
