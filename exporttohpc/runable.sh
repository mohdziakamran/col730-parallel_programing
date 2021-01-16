#!/bin/sh
echo "lets start"

g++ -fopenmp -c sort1.cpp sort2.cpp sort4.cpp sort6.cpp sort8.cpp sort10.cpp sort12.cpp sort14.cpp sort16.cpp sort18.cpp sort20.cpp 

g++ -fopenmp -o tester1 tester.cpp sort1.o
g++ -fopenmp -o tester2 tester.cpp sort2.o
g++ -fopenmp -o tester3 tester.cpp sort4.o
g++ -fopenmp -o tester4 tester.cpp sort6.o
g++ -fopenmp -o tester5 tester.cpp sort8.o
g++ -fopenmp -o tester6 tester.cpp sort10.o
g++ -fopenmp -o tester7 tester.cpp sort12.o
g++ -fopenmp -o tester8 tester.cpp sort14.o
g++ -fopenmp -o tester9 tester.cpp sort16.o
g++ -fopenmp -o tester10 tester.cpp sort18.o
g++ -fopenmp -o tester11 tester.cpp sort20.o

echo "for 1 thread"
./tester1
echo "for 2 threads"
./tester2
echo "for 4 threads"
./tester3
echo "for 6 threads"
./tester4
echo "for 8 threads"
./tester5
echo "for 10 threads"
./tester6
echo "for 12 threads"
./tester7
echo "for 14 threads"
./tester8
echo "for 16 threads"
./tester9
echo "for 18 threads"
./tester10
echo "for 20 threads"
./tester11

rm tester1 tester2 tester3 tester4 tester5 tester6 tester7 tester8 tester9 tester10

rm *.o