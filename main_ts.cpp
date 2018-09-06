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
    // argv[8] snapshot based sampling interval
    int snapshot_blocksize = atoi(argv[8]);
    // argv[9] interpolation option
    // 0 - TRICUBIC
    // else - TRILINEAR
    int option = atoi(argv[9]);
    if(option == 0) option = TRICUBIC;
    else option = TRILINEAR;
    // argv[10] mode
    // 0 - SZ snapshot & time
    // 1 - SZ snapshot & decimation time
    // 2 - decimation & SZ time
    // 3 - decimation snapshot & time
    int mode = atoi(argv[10]);
    size_t out_size = 0;

    switch(mode){
        case 0:
            SZ_compression_in_time<float>(argv[1], snapshot_num, interval, eb, n1, n2, n3, &out_size);
            SZ_decompression_in_time<float>(argv[1], snapshot_num, interval, n1, n2, n3);
            break;
        case 1:
            SZ_compress_snapshot_and_decimation_in_time<float>(argv[1], snapshot_num, interval, eb, n1, n2, n3);
            SZ_decompress_snapshot_and_interpolate_in_time<float>(argv[1], snapshot_num, interval, n1, n2, n3, option);
            break;
        case 2:
            decimation_snapshot_and_SZ_compression_in_time<float>(argv[1], snapshot_num, interval, snapshot_blocksize, eb, n1, n2, n3, option, &out_size);
            interpolate_snapshot_and_SZ_decompression_in_time<float>(argv[1], snapshot_num, interval, snapshot_blocksize, n1, n2, n3, option);
            break;
        case 3:
            decimation_sample_in_time_and_space<float>(argv[1], snapshot_num, interval, snapshot_blocksize, n1, n2, n3);
            decimation_interpolate_in_time_and_space<float>(argv[1], snapshot_num, interval, snapshot_blocksize, n1, n2, n3, option);
            break;
        default:
            std::cout << "No such mode!\n";
            exit(0);
    }

}