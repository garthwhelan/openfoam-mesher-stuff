main: build/main.o build/mesh.o 0 polyMesh build/NACA4.o
	g++ -o build/mesher build/main.o build/mesh.o build/NACA4.o

build/main.o: src/main.cpp build
	g++ -o build/main.o -c src/main.cpp

build/mesh.o: src/mesh.cpp build
	g++ -o build/mesh.o -c src/mesh.cpp

build/NACA4.o: src/mesh.cpp build
	g++ -o build/NACA4.o -c src/NACA4.cpp

0:
	mkdir 0

polyMesh:
	mkdir polyMesh

build:
	mkdir build

clean:
	rm -rf build polyMesh/* 0/*

build/test.o: src/test.cpp build
	g++ -o build/test.o -c src/test.cpp

test: build/test.o build/mesh.o
	g++ -o build/test build/test.o build/mesh.o -lboost_unit_test_framework && ./build/test
