import numpy as np
import os

def get_total_rate_distortion(nrmse, cr):
	total_nrmse = np.sqrt(np.mean(nrmse**2, axis=1))
	total_psnr = -20 * np.log10(total_nrmse)
	total_br = np.mean(32.0 / cr, axis=1)
	return total_br, total_psnr

def get_psnr_and_nrmse(data, dec_data):
    data_range = np.max(data) - np.min(data)
    diff = data - dec_data
    rmse = np.sqrt(np.mean(diff**2))
    if rmse == 0:
    	psnr = np.inf
    else:
	    psnr = 20 * np.log10(data_range / rmse)
    nrmse = rmse / data_range
    return psnr, nrmse

def get_statistics_1D(var, snapshot_num, interval, mode):
	interval_num = (snapshot_num - 1) // interval
	index = 1
	actual_snapshot = interval_num * interval
	br = np.zeros([actual_snapshot])
	nrmse = np.zeros([actual_snapshot])
	psnr = np.zeros([actual_snapshot])
	for i in range(actual_snapshot):
		data = np.fromfile("{}{:02d}.bin.dat".format(var, index), dtype=np.float32)[32:-64]
		dec_data = np.fromfile("{}{:02d}.bin.dat.{}.out".format(var, index, mode), dtype=np.float32)[32:-64]
		origin_size = os.path.getsize("{}{:02d}.bin.dat".format(var, index))
		compressed_size = os.path.getsize("{}{:02d}.bin.dat.{}".format(var, index, mode))
		br[i] = compressed_size * 32.0 / origin_size
		psnr[i], nrmse[i] = get_psnr_and_nrmse(data, dec_data)
		index += 1
	return br, psnr, nrmse

def get_statistics_3D(var, snapshot_num, interval, mode):
	interval_num = (snapshot_num - 1) // interval
	index = 1
	actual_snapshot = interval_num * interval
	br = np.zeros([actual_snapshot])
	nrmse = np.zeros([actual_snapshot])
	psnr = np.zeros([actual_snapshot])
	for i in range(actual_snapshot):
		data = np.fromfile("{}{:02d}.bin.dat".format(var, index), dtype=np.float32).reshape([100, 500, 500])[5:90, 5:490, 5:490]
		dec_data = np.fromfile("{}{:02d}.bin.dat.{}.out".format(var, index, mode), dtype=np.float32).reshape([100, 500, 500])[5:90, 5:490, 5:490]
		origin_size = os.path.getsize("{}{:02d}.bin.dat".format(var, index))
		compressed_size = os.path.getsize("{}{:02d}.bin.dat.{}".format(var, index, mode))
		br[i] = compressed_size * 32.0 / origin_size
		psnr[i], nrmse[i] = get_psnr_and_nrmse(data, dec_data)
		index += 1
	return cr, psnr, nrmse

def interpolation_in_space_1D(index, option, mode, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/EXAALT_multisteps/drbsd_test/exaalt"
	variables = np.array(["exaalt-x-", "exaalt-y-", "exaalt-z-"])
	modes = np.array(["szst", "szsdt", "dsszt", "dst"])
	var = variables[index]
	interval = 1
	interval_num = (83 - 1) // interval
	actual_snapshot = interval_num * interval
	sampling_dist = np.array([2, 4, 6, 8, 12, 16, 24, 32])
	cr = np.zeros([sampling_dist.size, actual_snapshot])
	psnr = np.zeros([sampling_dist.size, actual_snapshot])
	nrmse = np.zeros([sampling_dist.size, actual_snapshot])
	for i in range(sampling_dist.size):
		os.system("{} {}/{} 1 1 1077290 83 1 1e-3 {} {} {}".format(executable, directory, var, sampling_dist[i], option, mode))
		cr[i, :], psnr[i, :], nrmse[i, :] = get_statistics_1D("{}/{}".format(directory, var), 83, 1, modes[mode])
	np.savetxt("{}_{}_{}_snapshot_cr.txt".format(var, option, mode), cr)
	np.savetxt("{}_{}_{}_snapshot_psnr.txt".format(var, option, mode), psnr)
	np.savetxt("{}_{}_{}_snapshot_nrmse.txt".format(var, option, mode), nrmse)

def sz_space_1D(index, option=0, mode=0, executable="/home/xin/codes/test_compression/compression_ts"):
    directory = "/lcrc/project/ECP-EZ/public/compression/EXAALT_multisteps/drbsd_test/exaalt"
    variables = np.array(["exaalt-x-", "exaalt-y-", "exaalt-z-"])
    modes = np.array(["szst", "szsdt", "dsszt", "dst"])
    var = variables[index]
    interval = 1
    interval_num = (83 - 1) // interval
    actual_snapshot = interval_num * interval
    sampling_dist = np.array(["1e-1", "1e-2", "1e-3", "1e-4", "1e-5"])
    cr = np.zeros([sampling_dist.size, actual_snapshot])
    psnr = np.zeros([sampling_dist.size, actual_snapshot])
    nrmse = np.zeros([sampling_dist.size, actual_snapshot])
    for i in range(sampling_dist.size):
        os.system("{} {}/{} 1 1 1077290 83 1 {} 1 {} {}".format(executable, directory, var, sampling_dist[i], option, mode))
        cr[i, :], psnr[i, :], nrmse[i, :] = get_statistics_1D("{}/{}".format(directory, var), 83, 1, modes[mode])
    np.savetxt("{}_{}_{}_snapshot_cr.txt".format(var, option, mode), cr)
    np.savetxt("{}_{}_{}_snapshot_psnr.txt".format(var, option, mode), psnr)
    np.savetxt("{}_{}_{}_snapshot_nrmse.txt".format(var, option, mode), nrmse)

def interpolation_in_time_1D(index, interval, mode, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/EXAALT_multisteps/drbsd_test/exaalt"
	variables = np.array(["exaalt-x-", "exaalt-y-", "exaalt-z-"])
	modes = np.array(["szst", "szsdt", "dsszt", "dst"])
	var = variables[index]
	os.system("{} {}/{} 1 1 1077290 83 {} 0 1 0 {}".format(executable, directory, var, interval, mode))
	cr, psnr, nrmse = get_statistics_1D("{}/{}".format(directory, var), 83, interval, modes[mode])
	np.savetxt("{}_{}_time_interval_{}_cr.txt".format(var, modes[mode], interval), cr)
	np.savetxt("{}_{}_time_interval_{}_psnr.txt".format(var, modes[mode], interval), psnr)
	np.savetxt("{}_{}_time_interval_{}_nrmse.txt".format(var, modes[mode], interval), nrmse)

def sz_in_time_1D(index, eb, mode, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/EXAALT_multisteps/drbsd_test/exaalt"
	variables = np.array(["exaalt-x-", "exaalt-y-", "exaalt-z-"])
	modes = np.array(["szst", "szsdt", "dsszt", "dst"])
	var = variables[index]
	os.system("{} {}/{} 1 1 1077290 83 82 {} 1 0 {}".format(executable, directory, var, eb, mode))
	cr, psnr, nrmse = get_statistics_1D("{}/{}".format(directory, var), 83, 82, modes[mode])
	np.savetxt("{}_{}_time_eb_{}_cr.txt".format(var, modes[mode], eb), cr)
	np.savetxt("{}_{}_time_eb_{}_psnr.txt".format(var, modes[mode], eb), psnr)
	np.savetxt("{}_{}_time_eb_{}_nrmse.txt".format(var, modes[mode], eb), nrmse)

def interpolation_in_space_3D(index, option, mode, executable="/home/xin/codes/test_compression/compression_ts"):
	directory = "/lcrc/project/ECP-EZ/public/compression/test_data/Hurricane"
	variables = np.array(["QCLOUDf", "QGRAUPf", "QICEf", "QRAINf", "QSNOWf", "QVAPORf", "PRECIPf", "CLOUDf", "TCf", "Pf", "Uf", "Vf", "Wf"])
	modes = np.array(["szst", "szsdt", "dsszt", "dst"])
	var = variables[index]
	interval = 1
	interval_num = (83 - 1) // interval
	actual_snapshot = interval_num * interval
	sampling_dist = np.array([2, 3, 4, 5, 6])
	cr = np.zeros([sampling_dist.size, actual_snapshot])
	psnr = np.zeros([sampling_dist.size, actual_snapshot])
	nrmse = np.zeros([sampling_dist.size, actual_snapshot])
	total_br = np.zeros([sampling_dist.size])
	total_psnr = np.zeros([sampling_dist.size])
	for i in range(sampling_dist.size):
		os.system("{} {}/{} 100 500 500 48 1 0 {} {} {}".format(executable, directory, var, sampling_dist[i], option, mode))
		cr[i, :], psnr[i, :], nrmse[i, :] = get_statistics_3D("{}/{}".format(directory, var), 48, 1, modes[mode])
	np.savetxt("{}_{}_{}_snapshot_cr.txt".format(var, option, mode), cr)
	np.savetxt("{}_{}_{}_snapshot_psnr.txt".format(var, option, mode), psnr)
	np.savetxt("{}_{}_{}_snapshot_nrmse.txt".format(var, option, mode), nrmse)

