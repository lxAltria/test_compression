#ifndef _time_based_compression_h
#define _time_based_compression_h
#include "utils.h"
#include "sz.h"

template<typename Type>
unsigned char * SZ_compress_snapshot_based(Type * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size);

template<>
unsigned char * SZ_compress_snapshot_based<float>(float * ori_data, int n1, int n2, int n3, double error_bound, size_t * out_size){
	// int status = SZ_SCES;
	// float valueRangeSize = 0, medianValue = 0;
	// double realPrecision = getRealPrecision_float(valueRangeSize, REL, 0, error_bound, &status);
	// size_t dataLength = n1 * n2 * n3;
	// float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	// float max = min+valueRangeSize;
	// TightDataPointStorageF* tdps = SZ_compress_float_3D_MDQ(oriData, n1, n2, n3, realPrecision, valueRangeSize, medianValue);
	// unsigned char * comp_data;
	// convertTDPStoFlatBytes_float(tdps, &comp_data, out_size);
	// free_TightDataPointStorageF(tdps);
	unsigned char * comp_data = SZ_compress_args(SZ_FLOAT, ori_data, out_size, REL, error_bound, 0, 0, 0, 0, n1, n2, n3);
	return comp_data;
}

template<typename Type>
unsigned char * SZ_compress_time_based(Type * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size);

template<>
unsigned char * SZ_compress_time_based<float>(float * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size){
	int status = SZ_SCES;
	float valueRangeSize = 0, medianValue = 0;
	double realPrecision = getRealPrecision_float(valueRangeSize, REL, 0, error_bound, &status);
	size_t dataLength = n1 * n2 * n3;
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ_ts(oriData, dataLength, multisteps, realPrecision, valueRangeSize, medianValue);;
	unsigned char * comp_data;
	convertTDPStoFlatBytes_float(tdps, &comp_data, out_size);
	free_TightDataPointStorageF(tdps);
	return comp_data;
}


template<typename Type>
void SZ_compression_in_time(char * filename, int snapshot_num, int interval, double error_bound, int n1, int n2, int n3, size_t * out_size){
	SZ_Init(NULL);
	confparams_cpr->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_cpr->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	size_t dataLength = n1 * n2 * n3;
	Type * ori_data = (Type *) malloc(dataLength*sizeof(Type));

	int interval_num = (snapshot_num - 1)/ interval;
	TightDataPointStorageF* tdps = NULL;
	unsigned char * comp_data;
	size_t num_element;
	char filename_tmp[200];
	size_t total_size = 0;
	size_t index = 1;
	for(int i=0; i<interval_num; i++){
		size_t out_size;
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
		std::cout << "Interval " << i << ": " << filename_tmp <<  std::endl;
		// compress the first snapshot
		readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		comp_data = SZ_compress_snapshot_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
		writefile(strcat(filename_tmp, ".comp"), comp_data, out_size);
		total_size += out_size;
		free(comp_data);
		// compress the following interval-1 snapshot in time
		for(int j=0; j<interval-1; j++){
			std::cout << "snapshot " << j << ": " << filename_tmp <<  std::endl;
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
			readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
			comp_data = SZ_compress_time_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
			writefile(strcat(filename_tmp, ".comp"), comp_data, out_size);
			total_size += out_size;
			free(comp_data);
		}
	}
	SZ_Finalize();
	free(ori_data);
	*out_size = total_size;
}

// template<>
// void SZ_register_variable<float>(char * var_name, float ** data, double error_bound, int n1, int n2, int n3){
// 	SZ_registerVar(var_name, SZ_FLOAT, *data, REL, 0, error_bound, 0, 0, 0, n1, n2, n3);
// }

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