all: multiloop

CXXFLAGS=-Wall -g -pedantic -std=c++17 -I. -I/usr/include/eigen3

multiloop: multiloop.o lsq-plane.o harmonic.o triangle.o
	g++ -o $@ $^ -L. -lgeom
