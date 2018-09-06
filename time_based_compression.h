#ifndef _time_based_compression_h
#define _time_based_compression_h
#include "utils.h"
#include "sz.h"

template<typename Type>
unsigned char * SZ_compress_snapshot_based(Type * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size);

template<>
unsigned char * SZ_compress_snapshot_based<float>(float * ori_data, int n1, int n2, int n3, double error_bound, size_t * out_size){
	unsigned char * comp_data = SZ_compress_args(SZ_FLOAT, ori_data, out_size, REL, error_bound, error_bound, error_bound, 0, 0, n1, n2, n3);
	return comp_data;
}

template<typename Type>
unsigned char * SZ_compress_time_based(Type * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size);

template<>
unsigned char * SZ_compress_time_based<float>(float * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size){
	int status = SZ_SCES;
	float valueRangeSize = 0, medianValue = 0;
	size_t dataLength = n1 * n2 * n3;
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	double realPrecision = getRealPrecision_float(valueRangeSize, REL, 0, error_bound, &status);
	TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ_ts(oriData, dataLength, multisteps, realPrecision, valueRangeSize, medianValue);;
	unsigned char * comp_data;
	convertTDPStoFlatBytes_float(tdps, &comp_data, out_size);
	free_TightDataPointStorageF(tdps);
	if(confparams_cpr->szMode==SZ_BEST_SPEED){
		return comp_data;
	}
	else if(confparams_cpr->szMode==SZ_BEST_COMPRESSION || confparams_cpr->szMode==SZ_DEFAULT_COMPRESSION || confparams_cpr->szMode==SZ_TEMPORAL_COMPRESSION)
	{
		unsigned char * l_comp_data;
		*out_size = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, comp_data, *out_size, &l_comp_data);
		free(comp_data);
		return l_comp_data;
	}
	else
	{
		printf("Error: Wrong setting of confparams_cpr->szMode in the float compression.\n");
		status = SZ_MERR; //mode error			
	}
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
	if(sz_tsc == NULL) initSZ_TSC();
	multisteps = (sz_multisteps*)malloc(sizeof(sz_multisteps));
	memset(multisteps, 0, sizeof(sz_multisteps));
	multisteps->hist_data = (Type*)malloc(sizeof(Type)*dataLength);
	for(int i=0; i<interval_num; i++){
		size_t out_size;
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
		std::cout << "Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		// compress the first snapshot
		readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		comp_data = SZ_compress_snapshot_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
		writefile(strcat(filename_tmp, ".comp"), comp_data, out_size);
		total_size += out_size;
		free(comp_data);
		// compress the following interval-1 snapshot in time
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
			std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
			comp_data = SZ_compress_time_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
			writefile(strcat(filename_tmp, ".comp"), comp_data, out_size);
			total_size += out_size;
			free(comp_data);
		}
	}
	free_multisteps(multisteps);
	// if(sz_tsc) 
	SZ_Finalize();
	free(ori_data);
	*out_size = total_size;
}

template<typename Type>
void SZ_decompression_in_time(char * filename, int snapshot_num, int interval, int n1, int n2, int n3){
	// SZ_Init(NULL);
	if(confparams_dec==NULL)
		confparams_dec = (sz_params*)malloc(sizeof(sz_params));
	memset(confparams_dec, 0, sizeof(sz_params));

	confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_dec->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	size_t dataLength = n1 * n2 * n3;
	int interval_num = (snapshot_num - 1)/ interval;
	multisteps = (sz_multisteps*)malloc(sizeof(sz_multisteps));
	memset(multisteps, 0, sizeof(sz_multisteps));
	multisteps->hist_data = (float*)malloc(sizeof(Type)*dataLength);
	if(sz_tsc == NULL) initSZ_TSC();
	unsigned char * comp_data = (unsigned char *) malloc(sizeof(Type)*dataLength);
	// Type * dec_data = (Type *) malloc(sizeof(Type)*dataLength);
	size_t comp_data_size = 0;
	float * dec_data;
	float * ori_data = (float *) malloc(sizeof(float)*dataLength);
	char filename_tmp[200];
	size_t index = 1;
	for(int i=0; i<interval_num; i++){
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.comp", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat.comp", filename, index++);
		std::cout << "Decompression Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		readfile_to_buffer<unsigned char>(filename_tmp, &comp_data_size, comp_data);
		// decompress first snapshot
		dec_data = (float *)SZ_decompress(SZ_FLOAT, comp_data, comp_data_size, 0, 0, n1, n2, n3);
		memcpy(multisteps->hist_data, dec_data, sizeof(Type)*n1*n2*n3);
		writefile<float>(strcat(filename_tmp, ".out"), dec_data, n1*n2*n3);
		{
			// verify
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
			else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
			size_t num_element;
			readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
			verify(ori_data, dec_data, num_element, comp_data_size);
		}
		free(dec_data);
		confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.comp", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat.comp", filename, index++);
			std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			readfile_to_buffer<unsigned char>(filename_tmp, &comp_data_size, comp_data);
			TightDataPointStorageF* tdps;
			{
				int status = SZ_SCES;
				size_t dataLength = n1*n2*n3;
				
				//unsigned char* tmpBytes;
				size_t targetUncompressSize = dataLength <<2; //i.e., *4
				//tmpSize must be "much" smaller than dataLength
				size_t i, tmpSize = 8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
				unsigned char* szTmpBytes;	
				
				if(comp_data_size!=8+4+MetaDataByteLength && comp_data_size!=8+8+MetaDataByteLength) //4,8 means two posibilities of SZ_SIZE_TYPE
				{
					confparams_dec->losslessCompressor = is_lossless_compressed_data(comp_data, comp_data_size);
					if(confparams_dec->szMode!=SZ_TEMPORAL_COMPRESSION)
					{
						if(confparams_dec->losslessCompressor!=-1)
							confparams_dec->szMode = SZ_BEST_COMPRESSION;
						else
							confparams_dec->szMode = SZ_BEST_SPEED;			
					}
					
					if(confparams_dec->szMode==SZ_BEST_SPEED)
					{
						tmpSize = comp_data_size;
						szTmpBytes = comp_data;	
					}
					else if(confparams_dec->szMode==SZ_BEST_COMPRESSION || confparams_dec->szMode==SZ_DEFAULT_COMPRESSION || confparams_dec->szMode==SZ_TEMPORAL_COMPRESSION)
					{
						if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
							targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES; 
						tmpSize = sz_lossless_decompress(confparams_dec->losslessCompressor, comp_data, (unsigned long)comp_data_size, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE);//		(unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
						//szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
						//memcpy(szTmpBytes, tmpBytes, tmpSize);
						//free(tmpBytes); //release useless memory		
					}
					else
					{
						printf("Wrong value of confparams_dec->szMode in the double compressed bytes.\n");
						status = SZ_MERR;
						return;
					}	
				}
				else
					szTmpBytes = comp_data;
				//TODO: convert szTmpBytes to data array.
				int errBoundMode = new_TightDataPointStorageF_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
				if(confparams_dec->szMode!=SZ_BEST_SPEED && comp_data_size!=8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE)
					free(szTmpBytes);			
			}
			decompressDataSeries_float_1D_ts(&dec_data, n1*n2*n3, multisteps, tdps);	
			free_TightDataPointStorageF2(tdps);
			writefile<float>(strcat(filename_tmp, ".out"), dec_data, n1*n2*n3);
			{
				// verify
				if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
				else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
				size_t num_element;
				readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
				verify(ori_data, dec_data, num_element, comp_data_size);
			}
			free(dec_data);
		}
	}
	free_multisteps(multisteps);
	// SZ_Finalize();
	free(ori_data);
	free(comp_data);
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