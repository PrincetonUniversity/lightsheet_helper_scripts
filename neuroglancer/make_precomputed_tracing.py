#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:07:40 2019

@author: wanglab
"""

import os, numpy as np, igneous.task_creation as tc
from concurrent.futures import ProcessPoolExecutor
from PIL import Image
from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch
from taskqueue import LocalTaskQueue

shape = (813, 7166, 6046) # z,y,x

home_dir = '/home/wanglab/Documents/neuroglancer'

#tif_dir = '/home/wanglab/mounts/wang/pisano/tracing_output/antero_4x/20170115_tp_bl6_lob6b_ml_04/full_sizedatafld/20170115_tp_bl6_lob6b_ml_04_4x_647_008na_1hfds_z7d5um_75msec_10povlp_ch00'
tif_dir = '/home/wanglab/mounts/scratch/zmd/20170115_tp_bl6_lob6b_ml_04/cell_map_single_tifs'

def make_info_file(commit=True):
	info = CloudVolume.create_new_info(
		num_channels = 1,
		layer_type = 'image', # 'image' or 'segmentation'
		data_type = 'uint16', # 32 not necessary for Princeton atlas, but was for Allen atlas 
		encoding = 'raw', # other options: 'jpeg', 'compressed_segmentation' (req. uint32 or uint64)
		resolution = [ 1630, 1630, 7500 ], # X,Y,Z values in nanometers, 40 microns in each dim. 
		voxel_offset = [ 0, 0, 0 ], # values X,Y,Z values in voxels
		chunk_size = [ 1024, 1024, 1 ], # rechunk of image X,Y,Z in voxels, 
		volume_size = [6046,7166,813], # X,Y,Z size in voxels
	)

	# If you're using amazon or the local file system, you can replace 'gs' with 's3' or 'file'
	vol = CloudVolume('file:///home/wanglab/Documents/neuroglancer/20170115_tp_bl6_lob6b_ml_04/cells', info=info)
	vol.provenance.description = "TP NC timepoint"
	vol.provenance.owners = ['zmd@princeton.edu'] # list of contact email addresses
	if commit:
		vol.commit_info() # generates gs://bucket/dataset/layer/info json file
		vol.commit_provenance() # generates gs://bucket/dataset/layer/provenance json file
		print("Created CloudVolume info file: ",vol.info_cloudpath)
	return vol

def process(z):
	img_name = os.path.join(tif_dir, '20170115_tp_bl6_lob6b_ml_04_cell_map_Z%04d.tif' % z)
	print('Processing ', img_name)
	image = Image.open(img_name)
	width, height = image.size 
	array = np.array(list( image.getdata() ), dtype=np.uint16, order='F')
	array = array.reshape((1, height, width)).T
	vol[:,:, z] = array
	image.close()
	touch(os.path.join(progress_dir, str(z)))

def make_demo_downsample(mip_start=0,num_mips=3):
	cloudpath = 'file:///home/wanglab/Documents/neuroglancer/20170115_tp_bl6_lob6b_ml_04/cells'
	with LocalTaskQueue(parallel=8) as tq:
		tasks = tc.create_downsampling_tasks(
			cloudpath, 
			mip=mip_start, # Start downsampling from this mip level (writes to next level up)
			fill_missing=False, # Ignore missing chunks and fill them with black
			axis='z', 
			num_mips=num_mips, # number of downsamples to produce. Downloaded shape is chunk_size * 2^num_mip
			chunk_size=[ 128, 128, 32 ], # manually set chunk size of next scales, overrides preserve_chunk_size
			preserve_chunk_size=True, # use existing chunk size, don't halve to get more downsamples
		  )
		tq.insert_all(tasks)
	print("Done!")

#if __name__ == '__main__':
#
#	vol = make_info_file()
#	progress_dir = mkdir(home_dir + '/progress_20170115_tp_bl6_lob6b_ml_04_cells/') # unlike os.mkdir doesn't crash on prexisting 
#	done_files = set([ int(z) for z in os.listdir(progress_dir) ])
#	all_files = set(range(vol.bounds.minpt.z, vol.bounds.maxpt.z)) 
#	to_upload = [ int(z) for z in list(all_files.difference(done_files)) ]
#	to_upload.sort()
#
#
#	with ProcessPoolExecutor(max_workers=8) as executor:
#	    executor.map(process, to_upload)

if __name__ == '__main__':

    #make different resolutions
    make_demo_downsample(mip_start=0,num_mips=4)
