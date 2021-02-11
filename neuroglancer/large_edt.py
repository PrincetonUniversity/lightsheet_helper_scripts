import edt
import numpy as np


if __name__ == '__main__':
	labels = np.ones(shape=(687, 2560, 2160), dtype=np.uint16, order='F')
	dt = edt.edt(
	labels, anisotropy=(6, 6, 30), 
	black_border=True, order='F',
	parallel=4 # number of threads, <= 0 sets to num cpu
	) 
	print(labels)
