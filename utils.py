import numpy as np
import os

def get_psnr_and_nrmse(data, dec_data):
	data_range = np.max(data) - np.min(data)
	diff = data - dec_data
	rmse = np.sqrt(np.mean(diff**2))
	psnr = 20 * np.log10(data_range / rmse)
	nrmse = rmse / data_range
	return psnr, nrmse

def get_statistics(var, snapshot_num, interval, mode):
	interval_num = (snapshot_num - 1) // interval
	index = 1
	actual_snapshot = interval_num * interval
	cr = np.zeros([actual_snapshot])
	nrmse = np.zeros([actual_snapshot])
	psnr = np.zeros([actual_snapshot])
	for i in range(actual_snapshot):
		data = np.fromfile("{}{:02d}.bin.dat".format(var, index), dtype=np.float32)#.reshape([100, 500, 500])[10:80, 10:480, 10:480]
		dec_data = np.fromfile("{}{:02d}.bin.dat.{}.out".format(var, index, mode), dtype=np.float32)#.reshape([100, 500, 500])[10:80, 10:480, 10:480]
		origin_size = os.path.getsize("{}{:02d}.bin.dat".format(var, index))
		compressed_size = os.path.getsize("{}{:02d}.bin.dat.{}".format(var, index, mode))
		cr[i] = origin_size * 1.0 / compressed_size
		psnr[i], nrmse[i] = get_psnr_and_nrmse(data, dec_data)
		index += 1
	total_br = np.mean(32.0 / cr)
	total_nrmse = np.sqrt(np.mean(nrmse**2))
	total_psnr = -20 * np.log10(total_nrmse)
	return cr, psnr, nrmse, total_br, total_psnr

#./compression_ts /lcrc/project/ECP-EZ/public/compression/test_data/Hurricane/Uf 100 500 500 48 1 1e-3 4 0 0
def run_szst(interval, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/test_data/Hurricane"
	variables = np.array(["QCLOUDf", "QGRAUPf", "QICEf", "QRAINf", "QSNOWf", "QVAPORf", "PRECIPf", "CLOUDf", "TCf", "Pf", "Uf", "Vf", "Wf"])
	error_bounds = np.array(["1e-1", "1e-2", "1e-3", "1e-4"])
	interval_num = (48 - 1) // interval
	actual_snapshot = interval_num * interval
	mode = "szst"
	for i in range(variables.size):
		var = variables[i]
		cr = np.zeros([error_bounds.size, actual_snapshot])
		psnr = np.zeros([error_bounds.size, actual_snapshot])
		nrmse = np.zeros([error_bounds.size, actual_snapshot])
		total_br = np.zeros([error_bounds.size])
		total_psnr = np.zeros([error_bounds.size])
		for j in range(error_bounds.size):
			os.system("{} {}/{} 100 500 500 48 {} {} 0 0 0".format(executable, directory, var, interval, error_bounds[j]))
			cr[j, :], psnr[j, :], nrmse[j, :], total_br[j], total_psnr[j] = get_statistics("{}/{}".format(directory, var), 48, interval, mode)
		np.savetxt("{}_{}_{}_cr.txt".format(var, mode, interval), cr)
		np.savetxt("{}_{}_{}_psnr.txt".format(var, mode, interval), psnr)
		np.savetxt("{}_{}_{}_nrmse.txt".format(var, mode, interval), nrmse)

def run_szsdt(interval, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/test_data/Hurricane"
	variables = np.array(["QCLOUDf", "QGRAUPf", "QICEf", "QRAINf", "QSNOWf", "QVAPORf", "PRECIPf", "CLOUDf", "TCf", "Pf", "Uf", "Vf", "Wf"])
	error_bounds = np.array(["1e-1", "1e-2", "1e-3", "1e-4"])
	interval_num = (48 - 1) // interval
	actual_snapshot = interval_num * interval
	mode = "szsdt"
	for i in range(variables.size):
		var = variables[i]
		cr = np.zeros([error_bounds.size, actual_snapshot])
		psnr = np.zeros([error_bounds.size, actual_snapshot])
		nrmse = np.zeros([error_bounds.size, actual_snapshot])
		total_br = np.zeros([error_bounds.size])
		total_psnr = np.zeros([error_bounds.size])
		for j in range(error_bounds.size):
			os.system("{} {}/{} 100 500 500 48 {} {} 0 0 1".format(executable, directory, var, interval, error_bounds[j]))
			cr[j, :], psnr[j, :], nrmse[j, :], total_br[j], total_psnr[j] = get_statistics("{}/{}".format(directory, var), 48, interval, mode)
		np.savetxt("{}_{}_{}_cr.txt".format(var, mode, interval), cr)
		np.savetxt("{}_{}_{}_psnr.txt".format(var, mode, interval), psnr)
		np.savetxt("{}_{}_{}_nrmse.txt".format(var, mode, interval), nrmse)	

def run_dsszt(interval, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/test_data/Hurricane"
	variables = np.array(["QCLOUDf", "QGRAUPf", "QICEf", "QRAINf", "QSNOWf", "QVAPORf", "PRECIPf", "CLOUDf", "TCf", "Pf", "Uf", "Vf", "Wf"])
	error_bounds = np.array(["1e-1", "1e-2", "1e-3", "1e-4"])
	interval_in_space = np.array(["1", "2", "3", "4", "5", "6"])
	interval_num = (48 - 1) // interval
	actual_snapshot = interval_num * interval
	mode = "dsszt"
	for i in range(variables.size):
		var = variables[i]
		for p in range(2):
			for space in interval_in_space:
				cr = np.zeros([error_bounds.size, actual_snapshot])
				psnr = np.zeros([error_bounds.size, actual_snapshot])
				nrmse = np.zeros([error_bounds.size, actual_snapshot])
				total_br = np.zeros([error_bounds.size])
				total_psnr = np.zeros([error_bounds.size])
				for j in range(error_bounds.size):
					os.system("{} {}/{} 100 500 500 48 {} {} {} {} 2".format(executable, directory, var, interval, error_bounds[j], space, p))
					cr[j, :], psnr[j, :], nrmse[j, :], total_br[j], total_psnr[j] = get_statistics("{}/{}".format(directory, var), 48, interval, mode)
				np.savetxt("{}_{}_{}_{}_{}_cr.txt".format(var, mode, interval, space, p), cr)
				np.savetxt("{}_{}_{}_{}_{}_psnr.txt".format(var, mode, interval, space, p), psnr)
				np.savetxt("{}_{}_{}_{}_{}_nrmse.txt".format(var, mode, interval, space, p), nrmse)

def run_dst(interval, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/test_data/Hurricane"
	variables = np.array(["QCLOUDf", "QGRAUPf", "QICEf", "QRAINf", "QSNOWf", "QVAPORf", "PRECIPf", "CLOUDf", "TCf", "Pf", "Uf", "Vf", "Wf"])
	interval_in_space = np.array(["1", "2", "3", "4", "5", "6"])
	interval_num = (48 - 1) // interval
	actual_snapshot = interval_num * interval
	mode = "dst"
	for i in range(variables.size):
		var = variables[i]
		for p in range(2): 
			cr = np.zeros([interval_in_space.size, actual_snapshot])
			psnr = np.zeros([interval_in_space.size, actual_snapshot])
			nrmse = np.zeros([interval_in_space.size, actual_snapshot])
			total_br = np.zeros([interval_in_space.size])
			total_psnr = np.zeros([interval_in_space.size])
			for j in range(interval_in_space.size):
				os.system("{} {}/{} 100 500 500 48 {} 0 {} {} 3".format(executable, directory, var, interval, interval_in_space[j], p))
				cr[j, :], psnr[j, :], nrmse[j, :], total_br[j], total_psnr[j] = get_statistics("{}/{}".format(directory, var), 48, interval, mode)
			np.savetxt("{}_{}_{}_{}_cr.txt".format(var, mode, interval, p), cr)
			np.savetxt("{}_{}_{}_{}_psnr.txt".format(var, mode, interval, p), psnr)
			np.savetxt("{}_{}_{}_{}_nrmse.txt".format(var, mode, interval, p), nrmse)	

