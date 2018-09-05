#ifndef _sampling_h
#define _sampling_h

#include "utils.h"

template<typename Type>
Type * uniform_sampling(Type * data, int block_size, int n1, int n2, int n3, size_t * out_size);

#endif