#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:35:49 2019

@author: wanglab
"""


import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from PIL import Image

from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch

from taskqueue import LocalTaskQueue
import igneous.task_creation as tc

shape = (690, 2560, 2160) # z,y,x

home_dir = '/home/wanglab/Documents/neuroglancer'
# tif = '/jukebox/LightSheetTransfer/kelly/201908_cfos/190821_m61467_demons_20190702_1d3x_647_008na_1hfds_z5um_250msec_14-09-11/14-09-11_UltraII_raw_RawDataStack [00 x 00\]_C00_xyz-Table\ Z0156_UltraII\ Filter0000.ome.tif'
tif_dir = '/home/wanglab/mounts/scratch/zmd/an19/transformed_annotations/single_tifs'
# home_dir = '/home/ahoag/ngdemo'
# atlas_file = '/home/ahoag/mounts/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif'
# # processed_data_file = '/mounts/wang/Jess/lightsheet_output/201908_testing_ahoag/processed/an31/\
# an31_devcno_03082019_1d3x_488_017na_1hfds_z10um_100msec_resized_ch00_resampledforelastix.tif' # the result after step 2, with all z planes in one tif stack

def make_info_file(commit=True):
	info = CloudVolume.create_new_info(
		num_channels = 1,
		layer_type = 'image', # 'image' or 'segmentation'
		data_type = 'uint16', # 32 not necessary for Princeton atlas, but was for Allen atlas 
		encoding = 'raw', # other options: 'jpeg', 'compressed_segmentation' (req. uint32 or uint64)
		resolution = [ 5000, 5000, 1000 ], # X,Y,Z values in nanometers, 40 microns in each dim. 
		voxel_offset = [ 0, 0, 0 ], # values X,Y,Z values in voxels
		chunk_size = [ 1024, 1024, 1 ], # rechunk of image X,Y,Z in voxels, 
		volume_size = [2160,2560,690], # X,Y,Z size in voxels
	)

	# If you're using amazon or the local file system, you can replace 'gs' with 's3' or 'file'
	vol = CloudVolume('file:///home/wanglab/Documents/neuroglancer/an19/atlas', info=info)
	vol.provenance.description = "Jess DREADDs"
	vol.provenance.owners = ['zmd@princeton.edu'] # list of contact email addresses
	if commit:
		vol.commit_info() # generates gs://bucket/dataset/layer/info json file
		vol.commit_provenance() # generates gs://bucket/dataset/layer/provenance json file
		print("Created CloudVolume info file: ",vol.info_cloudpath)
	return vol

def process(z):
	img_name = os.path.join(tif_dir, 'an19_annotation_Z%04d.tif' % z)
	print('Processing ', img_name)
	image = Image.open(img_name)
	width, height = image.size 
	array = np.array(list( image.getdata() ), dtype=np.uint16, order='F')
	array = array.reshape((1, height, width)).T
	vol[:,:, z] = array
	image.close()
	touch(os.path.join(progress_dir, str(z)))

def make_demo_downsample(mip_start=0,num_mips=3):
	cloudpath = 'file:///home/wanglab/Documents/neuroglancer/an19/atlas'
	with LocalTaskQueue(parallel=8) as tq:
		tasks = tc.create_downsampling_tasks(
			cloudpath, 
			mip=mip_start, # Start downsampling from this mip level (writes to next level up)
			fill_missing=False, # Ignore missing chunks and fill them with black
			axis='z', 
			num_mips=num_mips, # number of downsamples to produce. Downloaded shape is chunk_size * 2^num_mip
			chunk_size=[ 128, 128, 8 ], # manually set chunk size of next scales, overrides preserve_chunk_size
			preserve_chunk_size=True, # use existing chunk size, don't halve to get more downsamples
		  )
		tq.insert_all(tasks)
	print("Done!")

if __name__ == '__main__':

#	vol = make_info_file()
#	progress_dir = mkdir(home_dir + '/progress_an19_atlas/') # unlike os.mkdir doesn't crash on prexisting 
#	done_files = set([ int(z) for z in os.listdir(progress_dir) ])
#	all_files = set(range(vol.bounds.minpt.z, vol.bounds.maxpt.z)) 
#	to_upload = [ int(z) for z in list(all_files.difference(done_files)) ]
#	to_upload.sort()
#
#
#	with ProcessPoolExecutor(max_workers=8) as executor:
#	    executor.map(process, to_upload)

#if __name__ == '__main__':

    #make different resolutions
    make_demo_downsample(mip_start=0,num_mips=3)
