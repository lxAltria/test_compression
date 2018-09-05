#include "compression.h"
#include "utils.h"
#include "time_based_compression.h"

int main(int argc, char ** argv){

    int status;
    size_t nbEle;
    // argv[1]: folder & var name
    // e.g. /lcrc/project/ECP-EZ/public/compression/Hurricane-ISABEL/clean-data-dtzip/Cloudf
    // argv[2-4] dimensions
    int n1 = atoi(argv[2]);
    int n2 = atoi(argv[3]);
    int n3 = atoi(argv[4]);
    // argv[5] snapshot_num
    int snapshot_num = atoi(argv[5]);
    // argv[6] interval
    int interval = atoi(argv[6]);
    // argv[7] eb in time/snapshot, may be separted later
    double eb = atof(argv[7]);
    size_t out_size = 0;
    SZ_compression_in_time<float>(argv[1], snapshot_num, interval, eb, n1, n2, n3, &out_size);
    SZ_decompression_in_time<float>(argv[1], snapshot_num, interval, n1, n2, n3);
    // verify(data, dec_data, nbEle, out_size);
    // free(comp_data);
    // free(data);
    // free(dec_data);
}