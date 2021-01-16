all:
	g++ -fopenmp -fPIC -c -Wall sort.cpp
	g++ -fopenmp -shared sort.o -o libpsort.so
clean:
	rm libpsort.so sort.o