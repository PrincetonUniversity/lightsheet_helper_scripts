import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from PIL import Image
import tifffile

from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch

from taskqueue import LocalTaskQueue
import igneous.task_creation as tc

home_dir = '/home/ahoag/ngdemo'
atlas_file = '/home/ahoag/mounts/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans.nrrd'

def make_info_file():
	info = CloudVolume.create_new_info(
		num_channels = 1,
		layer_type = 'segmentation', # 'image' or 'segmentation'
		data_type = 'uint32', # can pick any popular uint
		encoding = 'raw', # other options: 'jpeg', 'compressed_segmentation' (req. uint32 or uint64)
		resolution = [ 25000, 25000, 25000 ], # X,Y,Z values in nanometers, 40 microns in each dim
		voxel_offset = [ 0, 0, 0 ], # values X,Y,Z values in voxels
		chunk_size = [ 1024, 1024, 1 ], # rechunk of image X,Y,Z in voxels
		volume_size = [320, 528, 456], # X,Y,Z size in voxels
	)

	# If you're using amazon or the local file system, you can replace 'gs' with 's3' or 'file'
	vol = CloudVolume('file:///home/ahoag/ngdemo/demo_bucket/atlas/allenatlas_2017', info=info)
	vol.provenance.description = "Segmentation volume for the 3D labeled allen atlas"
	vol.provenance.owners = ['ahoag@princeton.edu'] # list of contact email addresses

	vol.commit_info() # generates gs://bucket/dataset/layer/info json file
	vol.commit_provenance() # generates gs://bucket/dataset/layer/provenance json file
	print("Created CloudVolume info file: ",vol.info_cloudpath)
	return vol


def process_slice(z):
	print('Processing slice z=',z)
	
	array = image[z-1].reshape((1,y_dim,x_dim)).T # the z-1 index is because the files in to_upload are 1-indexed

	vol[:,:, z] = array
	touch(os.path.join(progress_dir, str(z)))

def make_demo_mesh():
	# Mesh on 8 cores, use True to use all cores
	cloudpath = 'file:///home/ahoag/ngdemo/demo_bucket/atlas/allenatlas_2017'
	with LocalTaskQueue(parallel=8) as tq:
	  tasks = tc.create_meshing_tasks(cloudpath, mip=0, shape=(256, 256, 256))
	  tq.insert_all(tasks)
	  tasks = tc.create_mesh_manifest_tasks(cloudpath)
	  tq.insert_all(tasks)
	print("Done!")	

if __name__ == '__main__':
	""" Fill the CloudVolume() instance with data from the tif slices """
	vol = make_info_file()
	""" Now load the tifffile in its entirety """
	image = np.array(tifffile.imread(atlas_file),dtype=np.uint32, order='F') # F stands for fortran order
	z_dim,y_dim,x_dim = image.shape
	print(image.shape)

	progress_dir = mkdir(home_dir + '/progress_allenatlas_2017/') # unlike os.mkdir doesn't crash on prexisting 

	done_files = set([ int(z) for z in os.listdir(progress_dir) ])
	all_files = set(range(vol.bounds.minpt.z, vol.bounds.maxpt.z)) 

	to_upload = [ int(z) for z in list(all_files.difference(done_files)) ]
	to_upload.sort()
	print("Remaining slices to upload are:",to_upload)

	with ProcessPoolExecutor(max_workers=8) as executor:
	    executor.map(process_slice, to_upload)
	    



