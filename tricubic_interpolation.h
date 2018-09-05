#ifndef _tricubic_interpolation_h
#define _tricubic_interpolation_h

#include "utils.h"

template <typename Type>
void tricubic_interpolate_block(Type * _data, Type * _dec_data, int _spacing, int xi, int yi, int zi, int nx, int ny, int nz);

template <typename Type>
Type * tricubic_interpolation(Type * dec_comp_data, int block_size, int n3, int n2, int n1);

#endif