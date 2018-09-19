#ifndef _time_based_compression_h
#define _time_based_compression_h
#include "utils.h"
#include "compression.h"
#include "sz.h"
#include <ctime>
#include <sys/time.h>

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

void cost_start()
{
    totalCost = 0;
    gettimeofday(&costStart, NULL);
}

void cost_end()
{
    double elapsed;
    struct timeval costEnd;
    gettimeofday(&costEnd, NULL);
    elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
    totalCost += elapsed;
}

template<typename Type>
unsigned char * SZ_compress_snapshot_based(Type * oriData, int n1, int n2, int n3, double error_bound, size_t * out_size);

template<>
unsigned char * SZ_compress_snapshot_based<float>(float * ori_data, int n1, int n2, int n3, double error_bound, size_t * out_size){
	if(n1 == 1) n1 = 0;
	if(n2 == 1) n2 = 0;
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
	double elapsed_time = 0.0;
	for(int i=0; i<interval_num; i++){
		size_t out_size;
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
		// std::cout << "Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		// compress the first snapshot
		readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		cost_start();
		comp_data = SZ_compress_snapshot_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
		cost_end();
		elapsed_time += totalCost;
		writefile(strcat(filename_tmp, ".szst"), comp_data, out_size);
		total_size += out_size;
		free(comp_data);
		// compress the following interval-1 snapshot in time
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
			// std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
			cost_start();
			comp_data = SZ_compress_time_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
			cost_end();
			elapsed_time += totalCost;
			writefile(strcat(filename_tmp, ".szst"), comp_data, out_size);
			total_size += out_size;
			free(comp_data);
		}
	}
	free_multisteps(multisteps);
	// if(sz_tsc) 
	SZ_Finalize();
	free(ori_data);
	*out_size = total_size;
	std::cout << "SZST compression time: " << elapsed_time << " s, compression rate: " << dataLength * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
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
	if(n1 == 1) n1 = 0;
	if(n2 == 1) n2 = 0;
	int interval_num = (snapshot_num - 1)/ interval;
	multisteps = (sz_multisteps*)malloc(sizeof(sz_multisteps));
	memset(multisteps, 0, sizeof(sz_multisteps));
	multisteps->hist_data = (Type*)malloc(sizeof(Type)*dataLength);
	if(sz_tsc == NULL) initSZ_TSC();
	unsigned char * comp_data = (unsigned char *) malloc(sizeof(Type)*dataLength);
	// Type * dec_data = (Type *) malloc(sizeof(Type)*dataLength);
	size_t comp_data_size = 0;
	float * dec_data;
	Type * ori_data = (Type *) malloc(sizeof(Type)*dataLength);
	char filename_tmp[200];
	size_t index = 1;
	double elapsed_time = 0;
	for(int i=0; i<interval_num; i++){
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.szst", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat.szst", filename, index++);
		// std::cout << "Decompression Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		readfile_to_buffer<unsigned char>(filename_tmp, &comp_data_size, comp_data);
		// decompress first snapshot
		cost_start();
		dec_data = (float *)SZ_decompress(SZ_FLOAT, comp_data, comp_data_size, 0, 0, n1, n2, n3);
		memcpy(multisteps->hist_data, dec_data, sizeof(Type)*dataLength);
		cost_end();
		elapsed_time += totalCost;
		writefile<float>(strcat(filename_tmp, ".out"), dec_data, dataLength);
		// {
		// 	// verify
		// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
		// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
		// 	size_t num_element;
		// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		// 	verify(ori_data, dec_data, num_element, comp_data_size);
		// }
		free(dec_data);
		confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
		confparams_dec->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.szst", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat.szst", filename, index++);
			// std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			readfile_to_buffer<unsigned char>(filename_tmp, &comp_data_size, comp_data);
			cost_start();
			TightDataPointStorageF* tdps;
			{
				int status = SZ_SCES;				
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
				decompressDataSeries_float_1D_ts(&dec_data, dataLength, multisteps, tdps);	
				free_TightDataPointStorageF2(tdps);
				writefile<float>(strcat(filename_tmp, ".out"), dec_data, dataLength);
				if(confparams_dec->szMode!=SZ_BEST_SPEED && comp_data_size!=8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE)
					free(szTmpBytes);
				// {
				// 	// verify
				// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
				// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
				// 	size_t num_element;
				// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
				// 	verify(ori_data, dec_data, num_element, comp_data_size);
				// }
			}
			cost_end();
			elapsed_time += totalCost;
			free(dec_data);
		}
	}
	free_multisteps(multisteps);
	// SZ_Finalize();
	free(ori_data);
	free(comp_data);
	std::cout << "SZST decompression time: " << elapsed_time << " s, decompression rate: " << dataLength * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
}

// template<>
// void SZ_register_variable<float>(char * var_name, float ** data, double error_bound, int n1, int n2, int n3){
// 	SZ_registerVar(var_name, SZ_FLOAT, *data, REL, 0, error_bound, 0, 0, 0, n1, n2, n3);
// }

template<typename Type>
void decimation_sample_in_time_and_space(char * filename, int snapshot_num, int interval, int snapshot_blocksize, int n1, int n2, int n3){

	int interval_num = (snapshot_num - 1)/ interval;
	char filename_tmp[200];
	size_t index = 1;
	size_t num_element;
	Type * ori_data = (Type *) malloc(sizeof(Type)*n1*n2*n3);
	double elapsed_time = 0.0;
	for(int i=0; i<interval_num + 1; i++){
		size_t out_size;
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index);
		else sprintf(filename_tmp, "%s%d.bin.dat", filename, index);
		// std::cout << "Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		// compress the first snapshot
		readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		if(snapshot_blocksize > 1){
			cost_start();
			Type * comp_data = uniform_sampling(ori_data, snapshot_blocksize, n1, n2, n3, &out_size);
			cost_end();
			elapsed_time += totalCost;
			writefile(strcat(filename_tmp, ".dst"), comp_data, out_size/sizeof(Type));
			free(comp_data);
		}
		else{
			writefile(strcat(filename_tmp, ".dst"), ori_data, sizeof(Type)*num_element);
		}
		// skip interval_num snapshots
		index ++;
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
			writefile(strcat(filename_tmp, ".dst"), ori_data, 0);
			index ++;
		}
		// index += interval;
	}
	free(ori_data);
	std::cout << "DST compression time: " << elapsed_time << " s, compression rate: " << n1*n2*n3 * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
}

#define TRICUBIC 0
#define TRILINEAR 1
template<typename Type>
void decimation_interpolate_in_time_and_space(char * filename, int snapshot_num, int interval, int snapshot_blocksize, int n1, int n2, int n3, int option){
	int interval_num = (snapshot_num - 1)/ interval;
	char filename_tmp[200];
	size_t index = 1;
	size_t num_element;
	Type * comp_data = (Type *) malloc(sizeof(Type)*n1*n2*n3);
	Type ** dec_data = (Type **) malloc(sizeof(Type *)*(interval_num + 1));
	size_t comp_data_size;
	size_t total_size = 0;
	double elapsed_time = 0.0;
	for(int i=0; i<interval_num + 1; i++){
		// decompress the first snapshot
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.dst", filename, index);
		else sprintf(filename_tmp, "%s%d.bin.dat.dst", filename, index);
		readfile_to_buffer(filename_tmp, &comp_data_size, comp_data);
		if(snapshot_blocksize > 1){
			cost_start();
			if(option == TRICUBIC) dec_data[i] = tricubic_interpolation(comp_data, snapshot_blocksize, n1, n2, n3);
			else dec_data[i] = trilinear_interpolation(comp_data, snapshot_blocksize, n1, n2, n3);
			cost_end();
		}
		else{
			dec_data[i] = (Type *) malloc(sizeof(Type)*n1*n2*n3);
			memcpy(dec_data[i], comp_data, sizeof(Type)*n1*n2*n3);
		}
		elapsed_time += totalCost;
		index += interval;
		total_size += comp_data_size*sizeof(Type);
	}
	free(comp_data);
	index = 1;
	Type * dec_buffer = (Type *) malloc(sizeof(Type)*n1*n2*n3);
	Type * ori_data = (Type *) malloc(sizeof(Type)*n1*n2*n3);
	for(int i=0; i<interval_num; i++){
		// write data of sampled snapshot
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.dst", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat.dst", filename, index++);
		// std::cout << "Decompression Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		writefile<float>(strcat(filename_tmp, ".out"), dec_data[i], n1*n2*n3);
		// {
		// 	// verify
		// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
		// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
		// 	size_t num_element;
		// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		// 	verify(ori_data, dec_data[i], num_element, total_size/interval_num/interval);
		// }
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.dst", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat.dst", filename, index++);
			// std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			// linear interpolation
			cost_start();
			double dist = (j+1)*1.0/interval_num;
			for(int num=0; num<n1*n2*n3; num++){
				dec_buffer[num] = (1-dist)*dec_data[i][num] + dist*dec_data[i+1][num];
			}
			cost_end();
			elapsed_time += totalCost;
			writefile<float>(strcat(filename_tmp, ".out"), dec_buffer, n1*n2*n3);
			// {
			// 	// verify
			// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
			// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
			// 	size_t num_element;
			// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
			// 	verify(ori_data, dec_buffer, num_element, total_size/interval_num/interval);
			// }
		}
	}
	free(ori_data);
	free(dec_buffer);
	for(int i=0; i<interval_num + 1; i++){
		free(dec_data[i]);
	}
	free(dec_data);
	std::cout << "DST decompression time: " << elapsed_time << " s, decompression rate: " << n1*n2*n3 * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
}

template<typename Type>
void SZ_compress_snapshot_and_decimation_in_time(char * filename, int snapshot_num, int interval, double error_bound, int n1, int n2, int n3){

	SZ_Init(NULL);
	int interval_num = (snapshot_num - 1)/ interval;
	char filename_tmp[200];
	size_t index = 1;
	size_t num_element;
	Type * ori_data = (Type *) malloc(sizeof(Type)*n1*n2*n3);
	size_t total_size = 0;
	double elapsed_time = 0.0;
	for(int i=0; i<interval_num + 1; i++){
		size_t out_size;
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index);
		else sprintf(filename_tmp, "%s%d.bin.dat", filename, index);
		// std::cout << "Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		// compress the first snapshot
		readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		cost_start();
		unsigned char * comp_data = SZ_compress_snapshot_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
		cost_end();
		elapsed_time += totalCost;
		writefile(strcat(filename_tmp, ".szsdt"), comp_data, out_size);
		free(comp_data);
		// skip interval_num snapshots
		index ++;
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
			writefile(strcat(filename_tmp, ".szsdt"), ori_data, 0);
			index ++;
		}
		// index += interval;
		total_size += out_size;
	}
	free(ori_data);
	SZ_Finalize();
	std::cout << "SZSDT compression time: " << elapsed_time << " s, compression rate: " << n1*n2*n3 * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
}

template<typename Type>
void SZ_decompress_snapshot_and_interpolate_in_time(char * filename, int snapshot_num, int interval, int n1, int n2, int n3, int option){
	SZ_Init(NULL);
	int interval_num = (snapshot_num - 1)/ interval;
	char filename_tmp[200];
	size_t index = 1;
	size_t num_element;
	size_t dataLength = n1*n2*n3;
	if(n1 == 1) n1 = 0;
	if(n2 == 1) n2 = 0;
	unsigned char * comp_data = (unsigned char *) malloc(sizeof(Type)*dataLength);
	Type ** dec_data = (Type **) malloc(sizeof(Type *)*(interval_num + 1));
	size_t comp_data_size;
	size_t total_size = 0;
	double elapsed_time = 0;
	for(int i=0; i<interval_num + 1; i++){
		// decompress the first snapshot
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.szsdt", filename, index);
		else sprintf(filename_tmp, "%s%d.bin.dat.szsdt", filename, index);
		readfile_to_buffer(filename_tmp, &comp_data_size, comp_data);
		cost_start();
		dec_data[i] = (Type *)SZ_decompress(SZ_FLOAT, comp_data, comp_data_size, 0, 0, n1, n2, n3);
		cost_end();
		elapsed_time += totalCost;
		index += interval;
		total_size += comp_data_size;
	}
	free(comp_data);
	index = 1;
	Type * dec_buffer = (Type *) malloc(sizeof(Type)*dataLength);
	Type * ori_data = (Type *) malloc(sizeof(Type)*dataLength);
	for(int i=0; i<interval_num; i++){
		// write data of sampled snapshot
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.szsdt", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat.szsdt", filename, index++);
		// std::cout << "Decompression Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		writefile<float>(strcat(filename_tmp, ".out"), dec_data[i], dataLength);
		// {
		// 	// verify
		// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
		// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
		// 	size_t num_element;
		// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		// 	verify(ori_data, dec_data[i], num_element, total_size/interval_num/interval);
		// }
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.szsdt", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat.szsdt", filename, index++);
			// std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			// linear interpolation
			cost_start();
			double dist = (j+1)*1.0/interval_num;
			for(int num=0; num<dataLength; num++){
				dec_buffer[num] = (1-dist)*dec_data[i][num] + dist*dec_data[i+1][num];
			}
			cost_end();
			elapsed_time += totalCost;
			writefile<float>(strcat(filename_tmp, ".out"), dec_buffer, dataLength);
			// {
			// 	// verify
			// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
			// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
			// 	size_t num_element;
			// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
			// 	verify(ori_data, dec_buffer, num_element, total_size/interval_num/interval);
			// }
		}
	}
	free(ori_data);
	free(dec_buffer);
	for(int i=0; i<interval_num + 1; i++){
		free(dec_data[i]);
	}
	free(dec_data);
	SZ_Finalize();
	std::cout << "SZSDT decompression time: " << elapsed_time << " s, decompression rate: " << dataLength * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
}

template<typename Type>
void decimation_snapshot_and_SZ_compression_in_time(char * filename, int snapshot_num, int interval, int snapshot_blocksize, double error_bound, int n1, int n2, int n3, int option, size_t * out_size){
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
	double elapsed_time = 0.0;
	for(int i=0; i<interval_num; i++){
		size_t out_size;
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
		// std::cout << "Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		// compress the first snapshot
		readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		cost_start();
		comp_data = (unsigned char *) uniform_sampling(ori_data, snapshot_blocksize, n1, n2, n3, &out_size);
		cost_end();
		elapsed_time += totalCost;
		writefile(strcat(filename_tmp, ".dsszt"), comp_data, out_size);
		cost_start();
		{
			// decompression & record history
			Type * dec_data;
			if(option == TRICUBIC) dec_data = tricubic_interpolation((Type *)comp_data, snapshot_blocksize, n1, n2, n3);
			else dec_data = trilinear_interpolation((Type *)comp_data, snapshot_blocksize, n1, n2, n3);
			memcpy(multisteps->hist_data, dec_data, sizeof(Type)*dataLength);
			free(dec_data);
		}
		cost_end();
		elapsed_time += totalCost;
		total_size += out_size;
		free(comp_data);
		// compress the following interval-1 snapshot in time
		for(int j=0; j<interval-1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat", filename, index++);
			// std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
			cost_start();
			comp_data = SZ_compress_time_based<float>(ori_data, n1, n2, n3, error_bound, &out_size);
			cost_end();
			elapsed_time += totalCost;
			writefile(strcat(filename_tmp, ".dsszt"), comp_data, out_size);
			total_size += out_size;
			free(comp_data);
		}
	}
	free_multisteps(multisteps);
	// if(sz_tsc) 
	SZ_Finalize();
	free(ori_data);
	*out_size = total_size;
	std::cout << "DSSZT compression time: " << elapsed_time << " s, compression rate: " << n1*n2*n3 * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
}


template<typename Type>
void interpolate_snapshot_and_SZ_decompression_in_time(char * filename, int snapshot_num, int interval, int snapshot_blocksize, int n1, int n2, int n3, int option){
	// SZ_Init(NULL);
	if(confparams_dec==NULL)
		confparams_dec = (sz_params*)malloc(sizeof(sz_params));
	memset(confparams_dec, 0, sizeof(sz_params));
	if(exe_params==NULL)
		exe_params = (sz_exedata*)malloc(sizeof(sz_exedata));
	memset(exe_params, 0, sizeof(sz_exedata));

	confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
	confparams_dec->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	size_t dataLength = n1 * n2 * n3;
	int interval_num = (snapshot_num - 1)/ interval;
	multisteps = (sz_multisteps*)malloc(sizeof(sz_multisteps));
	memset(multisteps, 0, sizeof(sz_multisteps));
	multisteps->hist_data = (Type*)malloc(sizeof(Type)*dataLength);
	if(sz_tsc == NULL) initSZ_TSC();
	unsigned char * comp_data = (unsigned char *) malloc(sizeof(Type)*dataLength);
	// Type * dec_data = (Type *) malloc(sizeof(Type)*dataLength);
	size_t comp_data_size = 0;
	Type * dec_data;
	Type * ori_data = (Type *) malloc(sizeof(Type)*dataLength);
	char filename_tmp[200];
	size_t index = 1;
	size_t total_size = 0;
	double elapsed_time = 0.0;
	for(int i=0; i<interval_num; i++){
		if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.dsszt", filename, index++);
		else sprintf(filename_tmp, "%s%d.bin.dat.dsszt", filename, index++);
		// std::cout << "Decompression Interval " << i << ":\nsnapshot 0: " << filename_tmp <<  std::endl;
		readfile_to_buffer<unsigned char>(filename_tmp, &comp_data_size, comp_data);
		total_size += comp_data_size;
		cost_start();
		// decompress first snapshot
		{
			// decompression & record history
			if(option == TRICUBIC) dec_data = tricubic_interpolation((Type *)comp_data, snapshot_blocksize, n1, n2, n3);
			else dec_data = trilinear_interpolation((Type *)comp_data, snapshot_blocksize, n1, n2, n3);
			memcpy(multisteps->hist_data, dec_data, sizeof(Type)*dataLength);
		}
		cost_end();
		elapsed_time += totalCost;
		writefile<float>(strcat(filename_tmp, ".out"), dec_data, n1*n2*n3);
		// {
		// 	// verify
		// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
		// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
		// 	size_t num_element;
		// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
		// 	verify(ori_data, dec_data, num_element, comp_data_size);
		// }
		free(dec_data);
		confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
		confparams_dec->predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
		for(int j=0; j<interval - 1; j++){
			if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat.dsszt", filename, index++);
			else sprintf(filename_tmp, "%s%d.bin.dat.dsszt", filename, index++);
			// std::cout << "snapshot " << j+1 << ": " << filename_tmp <<  std::endl;
			readfile_to_buffer<unsigned char>(filename_tmp, &comp_data_size, comp_data);
			total_size += comp_data_size;
			cost_start();
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
				decompressDataSeries_float_1D_ts(&dec_data, n1*n2*n3, multisteps, tdps);	
				free_TightDataPointStorageF2(tdps);
				writefile<float>(strcat(filename_tmp, ".out"), dec_data, n1*n2*n3);
				if(confparams_dec->szMode!=SZ_BEST_SPEED && comp_data_size!=8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE)
					free(szTmpBytes);
				// {
				// 	// verify
				// 	if(index < 10) sprintf(filename_tmp, "%s0%d.bin.dat", filename, index - 1);
				// 	else sprintf(filename_tmp, "%s%d.bin.dat", filename, index - 1);
				// 	size_t num_element;
				// 	readfile_to_buffer<float>(filename_tmp, &num_element, ori_data);
				// 	verify(ori_data, dec_data, num_element, comp_data_size);
				// }
			}
			cost_end();
			elapsed_time += totalCost;
			free(dec_data);
		}
	}
	// std::cout << "Total compression ratio: " << n1*n2*n3*sizeof(Type) * 1.0 * (interval_num * interval) / total_size << std::endl; 
	free_multisteps(multisteps);
	// SZ_Finalize();
	free(ori_data);
	free(comp_data);
	std::cout << "DSSZT decompression time: " << elapsed_time << " s, decompression rate: " << n1*n2*n3 * sizeof(Type) * 1.0 * interval_num * interval / elapsed_time / 1024 / 1024 << " MB/s" << std::endl;
}

#endif
