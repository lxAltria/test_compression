import numpy as np

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

def interpolation_in_space_1D(index):
	directory = "/lcrc/project/ECP-EZ/public/compression/EXAALT_multisteps/drbsd_test/exaalt"
	variables = np.array(["exaalt-x-", "exaalt-y-", "exaalt-z-"])
	modes = np.array(["szst", "szsdt", "dsszt", "dst"])
	var = variables[index]
	interval_num = (83 - 1) // interval
	actual_snapshot = interval_num * interval
	sampling_dist = np.array([2, 4, 6, 8, 12, 16, 24, 32])
	cr = np.zeros([sampling_dist.size, actual_snapshot])
	psnr = np.zeros([sampling_dist.size, actual_snapshot])
	nrmse = np.zeros([sampling_dist.size, actual_snapshot])
	total_br = np.zeros([sampling_dist.size])
	total_psnr = np.zeros([sampling_dist.size])
	for i in range(sampling_dist.size):
		os.system("{} {}/{} 1 1 1077290 83 1 1e-3 {} {} {}".format(executable, directory, var, sampling_dist[j], option, mode))
		cr[i, :], psnr[i, :], nrmse[i, :], total_br[i], total_psnr[i] = get_statistics("{}/{}".format(directory, var), actual_snapshot, 1, modes[mode])
	np.savetxt("{}_{}_{}_snapshot_cr.txt".format(var, option, mode), cr)
	np.savetxt("{}_{}_{}_snapshot_psnr.txt".format(var, option, mode), psnr)
	np.savetxt("{}_{}_{}_snapshot_nrmse.txt".format(var, option, mode), nrmse)
