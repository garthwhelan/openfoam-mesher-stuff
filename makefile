main: main.o mesh.o
	g++ -o main main.o mesh.o

main.o: main.cpp
	g++ -c main.cpp

mesh.o: mesh.cpp
	g++ -c mesh.cpp

clean:
	rm -rf main main.o mesh.o polyMesh/*
