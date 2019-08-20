
main: build/main.o build/mesh.o
	g++ -o build/mesher build/main.o build/mesh.o

build/main.o: main.cpp build
	g++ -o build/main.o -c main.cpp

build/mesh.o: mesh.cpp build
	g++ -o build/mesh.o -c mesh.cpp

build:
	mkdir build

clean:
	rm -rf build polyMesh/*

build/test.o: test.cpp build
	g++ -o build/test.o -c test.cpp

test: test.o mesh.o
	g++ -o build/test build/test.o build/mesh.o -lboost_unit_test_framework && ./build/test
