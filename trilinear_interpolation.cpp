#include "trilinear_interpolation.h"

static inline int _get_index(int i, int j, int k, int n1, int n2, int n3, size_t strip0, size_t strip1){
    if((i %= n1) < 0) i += n1;
    if((j %= n2) < 0) j += n2;
    if((k %= n3) < 0) k += n3;
    return i*strip0 + j*strip1 + k;
}

static inline int _get_index_1D(int i, int n){
	if((i %= n) < 0) i += n;
	return i;
}

template<typename Type>
Type * trilinear_interpolation(Type * comp_data, int block_size, int n1, int n2, int n3){

	if(n1 > 1){
		int nx = (n1 - 1)/block_size + 1;
		int ny = (n2 - 1)/block_size + 1;
		int nz = (n3 - 1)/block_size + 1;
		int comp_strip0 = ny * nz;
		int comp_strip1 = nz;
		int dec_strip0 = n2*n3;
		int dec_strip1 = n3;
		Type * dec_data = (Type *) malloc(n1*n2*n3*sizeof(Type));
		Type * dec_data_pos = NULL;
		double c00, c01, c10, c11, c0, c1;
		double * coeff = (double *) malloc(block_size * sizeof(double));
		for(int i=0; i<block_size; i++){
			coeff[i] = i * 1.0 / block_size;
		}
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				for(int k=0; k<nz; k++){
	                int bsi = block_size, bsj = block_size, bsk = block_size;
	                if(i*block_size + block_size > n1) bsi = n1 - i*block_size - block_size;
	                if(j*block_size + block_size > n2) bsj = n2 - j*block_size - block_size;
	                if(k*block_size + block_size > n3) bsk = n3 - k*block_size - block_size;
	                dec_data_pos = dec_data + i*block_size*dec_strip0 + j*block_size*dec_strip1 + k*block_size;
					for(int ii=0; ii<bsi; ii++){
						for(int jj=0; jj<bsj; jj++){
							for(int kk=0; kk<bsk; kk++){
								c00 = comp_data[_get_index(i, j, k, nx, ny, nz, comp_strip0, comp_strip1)] * (1 - coeff[ii]) + comp_data[_get_index(i, j, k+1, nx, ny, nz, comp_strip0, comp_strip1)] * coeff[ii];
								c01 = comp_data[_get_index(i, j+1, k, nx, ny, nz, comp_strip0, comp_strip1)] * (1 - coeff[ii]) + comp_data[_get_index(i, j+1, k+1, nx, ny, nz, comp_strip0, comp_strip1)] * coeff[ii];
								c10 = comp_data[_get_index(i+1, j, k, nx, ny, nz, comp_strip0, comp_strip1)] * (1 - coeff[ii]) + comp_data[_get_index(i+1, j, k+1, nx, ny, nz, comp_strip0, comp_strip1)] * coeff[ii];
								c11 = comp_data[_get_index(i+1, j+1, k, nx, ny, nz, comp_strip0, comp_strip1)] * (1 - coeff[ii]) + comp_data[_get_index(i+1, j+1, k+1, nx, ny, nz, comp_strip0, comp_strip1)] * coeff[ii];
								c0 = c00 * (1 - coeff[jj]) + c10 * coeff[jj];
								c1 = c01 * (1 - coeff[jj]) + c11 * coeff[jj];
	                            *(dec_data_pos++) = c0 * (1 - coeff[ii]) + c1 * coeff[ii];
							}
	                        dec_data_pos += dec_strip1 - bsk;
						}
	                    dec_data_pos += dec_strip0 - bsj * dec_strip1;
					}				 
				}
			}
		}
		free(coeff);
		return dec_data;
	}
	else if(n1 == 1 && n2 == 1){
		int nz = (n3 - 1)/block_size + 1;
		Type * dec_data = (Type *) malloc(n3*sizeof(Type));
		Type * dec_data_pos = dec_data;
		double * coeff = (double *) malloc(block_size * sizeof(double));
		for(int i=0; i<block_size; i++){
			coeff[i] = i * 1.0 / block_size;
		}
		for(int i=0; i<nz; i++){
			int bsi = block_size;
            if(i*block_size + block_size > n3) bsi = n3 - i*block_size - block_size;
			for(int ii=0; ii<bsi; ii++){
				*(dec_data_pos++) = (1 - coeff[ii]) * comp_data[_get_index_1D(i, nz)] + coeff[ii] * comp_data[_get_index_1D(i+1, nz)];
			}
		}
		return dec_data;
	}
	else{
		printf("Not support dimensions other than 1 and 3\n");
		exit(0);
	}
}

template float * trilinear_interpolation<float>(float * data, int block_size, int n1, int n2, int n3);
