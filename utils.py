import numpy as np

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
		data = np.fromfile("{}{:02d}.bin.dat".format(var, index), dtype=np.float32)
		dec_data = np.fromfile("{}{:02d}.bin.dat.{}.out".format(var, index, mode), dtype=np.float32)
		origin_size = os.path.getsize("{}{:02d}.bin.dat".format(var, index))
		compressed_size = os.path.getsize("{}{:02d}.bin.dat.{}".format(var, index, mode))
		cr[i] = origin_size * 1.0 / compressed_size
		psnr[i], nrmse[i] = get_psnr_and_nrmse(data, dec_data)
	total_br = np.mean(32.0 / cr)
	total_nrmse = np.sqrt(np.mean(nrmse**2))
	total_psnr = -20 * np.log10(total_nrmse)
	return cr, psnr, nrmse, total_br, total_psnr

