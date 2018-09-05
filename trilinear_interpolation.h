#ifndef _trilinear_interpolation_h
#define _trilinear_interpolation_h

#include "utils.h"

template<typename Type>
Type * trilinear_interpolation(Type * comp_data, int block_size, int n1, int n2, int n3);

#endif