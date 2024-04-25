sph:
	gcc -c -Wall -Werror -fpic sph.c -lm -lgsl -o sph.o
	gcc -shared -o libsph.so sph.o -lgsl -lgslcblas -lm