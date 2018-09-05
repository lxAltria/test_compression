#ifndef _time_based_compression_h
#define _time_based_compression_h
#include "utils.h"
#include "sz.h"

template<typename Type>
void SZ_compression_in_time(char * filename, int snapshot_num, int interval, double error_bound, int n1, int n2, int n3){
	SZ_Init(NULL);
	confparams_cpr->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_cpr->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	size_t dataLength = n1 * n2 * n3;
	Type * data = (Type*)malloc(sizeof(Type)*dataLength);

	int interval_num = (snapshot_num - 1)/ interval;
	TightDataPointStorageF* tdps = NULL;
	unsigned char * comp_data;
	size_t out_size;
	size_t num_element;
	char filename_tmp[200];
	for(int i=0; i<interval_num; i++){
		sprintf("%s_step_%d.dat", filename, i);
		// compress the first snapshot
		oridata = readfile<float>(filename_tmp, &num_element);
		comp_data = SZ_compress_snapshot_based(oriData, n1, n2, n3, error_bound, &out_size);
		// compress the following interval-1 snapshot in time
		for(int j=0; j<interval-1; j++){
			comp_data = SZ_compress_time_based(oriData, n1, n2, n3, error_bound, &out_size);
		}
		writefile(strcat(filename_tmp, ".comp"), comp_data, out_size);
	}
	SZ_Finalize();
}

// template<>
// void SZ_register_variable<float>(char * var_name, float ** data, double error_bound, int n1, int n2, int n3){
// 	SZ_registerVar(var_name, SZ_FLOAT, *data, REL, 0, error_bound, 0, 0, 0, n1, n2, n3);
// }

template<>
unsigned char * SZ_compress_snapshot_based<float>(float * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size){
	int status = SZ_SCES;
	float valueRangeSize = 0, medianValue = 0;
	realPrecision = getRealPrecision_float(valueRangeSize, REL, 0, error_bound, &status);
	size_t dataLength = n1 * n2 * n3;
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	TightDataPointStorageF* tdps = SZ_compress_float_3D_MDQ(oriData, n1, n2, n3, realPrecision, valueRangeSize, medianValue_f);
	unsigned char * comp_data;
	convertTDPStoFlatBytes_float(tdps, &comp_data, out_size);
	free_TightDataPointStorageF(tdps);
	return return comp_data;
}

template<>
unsigned char * SZ_compress_time_based<float>(float * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size){
	int status = SZ_SCES;
	float valueRangeSize = 0, medianValue = 0;
	realPrecision = getRealPrecision_float(valueRangeSize, REL, 0, error_bound, &status);
	size_t dataLength = n1 * n2 * n3;
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ_ts(oriData, dataLength, multisteps, realPrecision, valueRangeSize, medianValue_f);;
	unsigned char * comp_data;
	convertTDPStoFlatBytes_float(tdps, &comp_data, out_size);
	free_TightDataPointStorageF(tdps);
	return return comp_data;
}

template<typename Type>
void decimation_sample_in_time(char * filename, int snapshot_num, double error_bound, int n1, int n2, int n3){

	// compress the first snapshot

	// done

}

template<typename Type>
void decimation_interpolate_in_time(char * filename, int snapshot_num, double error_bound, int n1, int n2, int n3){

	// decompress the first snapshot

	// interpolation in time

}

#endif