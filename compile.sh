#! /bin/bash

g++ -c tricubic_interpolation.cpp -o tci.o
g++ -c trilinear_interpolation.cpp -o tli.o
g++ -c sampling.cpp -o sample.o
g++ main.cpp sample.o tci.o tli.o -o compression
g++ main_ts.cpp sample.o tci.o tli.o -o compression_ts -I/Users/LiangXin/utils/sz_ts/include -L/Users/LiangXin/utils/sz_ts/lib -lSZ -lzlib -lzstd -lm

