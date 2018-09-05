#include "sampling.h"

template<typename Type>
Type * uniform_sampling(Type * data, int block_size, int n1, int n2, int n3, size_t * out_size){
    int nx = (n1 - 1)/block_size + 1;
    int ny = (n2 - 1)/block_size + 1;
    int nz = (n3 - 1)/block_size + 1;
    Type * comp_data = (Type *) malloc(nx*ny*nz*sizeof(Type));
    Type * comp_data_pos = comp_data;
    int strip_0 = n2*n3;
    int strip_1 = n3;
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                *(comp_data_pos++) = *(data + i*block_size*strip_0 + j*block_size*strip_1 + k*block_size);
            }
        }
    }
    *out_size = nx*ny*nz*sizeof(Type);
    return comp_data;
}

template float * uniform_sampling<float>(float * data, int block_size, int n1, int n2, int n3, size_t * out_size);
