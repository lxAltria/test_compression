#! /bin/bash

g++ -c tricubic_interpolation.cpp -o tci.o
g++ -c trilinear_interpolation.cpp -o tli.o
g++ -c sampling.cpp -o sample.o
g++ main.cpp sample.o tci.o tli.o -o compression

