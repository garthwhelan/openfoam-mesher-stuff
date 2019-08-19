
main: main.o mesh.o
	g++ -o main main.o mesh.o

main.o: main.cpp
	g++ -c main.cpp

mesh.o: mesh.cpp
	g++ -c mesh.cpp

clean:
	rm -rf main main.o mesh.o test.o test polyMesh/*

test.o: test.cpp
	g++ -c test.cpp

test: test.o mesh.o
	g++ -o test test.o mesh.o -lboost_unit_test_framework
