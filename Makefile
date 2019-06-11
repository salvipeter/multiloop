all: multiloop

CXXFLAGS=-Wall -g -pedantic -std=c++17 -I.

multiloop: multiloop.o harmonic.o triangle.o
	g++ -o $@ $^ -L. -lgeom
