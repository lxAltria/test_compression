#include "compression.h"
#include "utils.h"

int main(int argc, char ** argv){

    int status;
    size_t nbEle;
    float * data = readfile<float>(argv[1], &nbEle);
    int n1 = atoi(argv[2]);
    int n2 = atoi(argv[3]);
    int n3 = atoi(argv[4]);
    int bs = atoi(argv[5]);
    size_t out_size = 0;
    float * comp_data = uniform_sampling(data, bs, n1, n2, n3, &out_size);
    std::cout << "compression done\n";
    float * dec_data = trilinear_interpolation(comp_data, bs, n1, n2, n3);
    // float * dec_data = tricubic_interpolation(comp_data, bs, n1, n2, n3);
    std::cout << "decompression done\n";
    // writefile<float>(strcat(argv[1], ".out"), data, nbEle);
    verify(data, dec_data, nbEle, out_size);
    free(comp_data);
    free(data);
    free(dec_data);
}