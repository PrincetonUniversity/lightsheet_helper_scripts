#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:07:40 2019

@author: wanglab
"""

import os, numpy as np, igneous.task_creation as tc, tifffile, sys
from concurrent.futures import ProcessPoolExecutor
from PIL import Image
from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch
from taskqueue import LocalTaskQueue

def make_info_file(brain, home_dir, volume_size, type_vol = "647", commit=True):
	info = CloudVolume.create_new_info(
		num_channels = 1,
		layer_type = "image", # "image" or "segmentation"
		data_type = "uint16", # 32 not necessary for Princeton atlas, but was for Allen atlas 
		encoding = "raw", # other options: "jpeg", "compressed_segmentation" (req. uint32 or uint64)
		resolution = [ 1630, 1630, 7500 ], # X,Y,Z values in nanometers, 40 microns in each dim. 
		voxel_offset = [ 0, 0, 0 ], # values X,Y,Z values in voxels
		chunk_size = [ 1024, 1024, 1 ], # rechunk of image X,Y,Z in voxels, 
		volume_size = volume_size, # X,Y,Z size in voxels
	)

	# If you"re using amazon or the local file system, you can replace "gs" with "s3" or "file"
	vol = CloudVolume("file://"+home_dir+"/"+brain+"/"+type_vol, info=info)
	vol.provenance.description = "TP tracing"
	vol.provenance.owners = ["zmd@princeton.edu"] # list of contact email addresses
	if commit:
		vol.commit_info() # generates gs://bucket/dataset/layer/info json file
		vol.commit_provenance() # generates gs://bucket/dataset/layer/provenance json file
		print("Created CloudVolume info file: ",vol.info_cloudpath)
	return vol

def process(z):
    
    if "full_sizedatafld" in tif_dir: #if raw images
        img_name = os.path.join(tif_dir, os.path.basename(tif_dir)+"_C00_Z%04d.tif" % z)
    elif "transformed_annotations" in tif_dir: #if atlas vol
        brain = os.path.basename(os.path.dirname(os.path.dirname(tif_dir)))
        img_name = os.path.join(tif_dir, brain+"_annotation_Z%04d.tif" % z)
    else:
        img_name = os.path.join(tif_dir, "cells_%04d.tif" % z)
        
    print("Processing ", img_name)
    image = Image.open(img_name)
    width, height = image.size
    array = np.array(list( image.getdata() ), dtype=np.uint16, order="F")
    array = array.reshape((1, height, width)).T
    vol[:,:, z] = array
    image.close()
    touch(os.path.join(progress_dir, str(z)))

def make_demo_downsample(type_vol="647", mip_start=0, num_mips=3):
	cloudpath = "file://"+home_dir+"/"+brain+"/"+type_vol
	with LocalTaskQueue(parallel=8) as tq:
		tasks = tc.create_downsampling_tasks(
			cloudpath, 
			mip=mip_start, # Start downsampling from this mip level (writes to next level up)
			fill_missing=False, # Ignore missing chunks and fill them with black
			axis="z", 
			num_mips=num_mips, # number of downsamples to produce. Downloaded shape is chunk_size * 2^num_mip
			chunk_size=[ 128, 128, 32 ], # manually set chunk size of next scales, overrides preserve_chunk_size
			preserve_chunk_size=True, # use existing chunk size, don"t halve to get more downsamples
		  )
		tq.insert_all(tasks)
	print("Done!")

if __name__ == "__main__":
    
    #setting dirs
    home_dir = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125"
    
    brain = str(sys.argv[1])
    type_vol = str(sys.argv[2]) #to separate raw images from atlas
    
    if type_vol == "647":
        src = "/jukebox/wang/pisano/tracing_output/antero_4x"
        fullszdt_dir = os.path.join(os.path.join(src, brain), "full_sizedatafld")
        tif_dir = [os.path.join(fullszdt_dir, xx) for xx in os.listdir(fullszdt_dir) if "647" in xx][0]
    elif type_vol == "atlas":
        src = "/jukebox/scratch/zmd"
        tif_dir = os.path.join(os.path.join(src, brain), "transformed_annotations/single_tifs")
    elif type_vol == "cells":
        src = "/jukebox/scratch/zmd"
        tif_dir = os.path.join(os.path.join(src, brain), "cells")
        
    #get x,y,z resolution
    y,x = tifffile.imread(os.path.join(tif_dir, os.listdir(tif_dir)[0])).shape
    volume_size = [x, y, len(os.listdir(tif_dir))]
    vol = make_info_file(brain, home_dir, volume_size, type_vol = type_vol)
    
    progress_dir = mkdir(home_dir + "/progress_"+brain+"_"+type_vol)
    done_files = set([ int(z) for z in os.listdir(progress_dir) ])
    all_files = set(range(vol.bounds.minpt.z, vol.bounds.maxpt.z))
    to_upload = [ int(z) for z in list(all_files.difference(done_files)) ]
    to_upload.sort()
    
    with ProcessPoolExecutor(max_workers=12) as executor:
        executor.map(process, to_upload)
    
    make_demo_downsample(type_vol, mip_start=0,num_mips=5)