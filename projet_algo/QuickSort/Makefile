all: quickSort

quickSort: quick_sort.o
	gcc -Wall -o quickSort quick_sort.o

main.o: quick_sort.c quick_sort.h
	gcc -Wall -c -o quick_sort.o quick_sort.c

clean:
	rm *.o *.out

mrproper: clean
	rm -rf main.out
