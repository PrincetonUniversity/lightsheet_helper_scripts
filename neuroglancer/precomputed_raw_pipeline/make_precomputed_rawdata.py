#! /bin/env python

import os, sys
import glob
from concurrent.futures import ProcessPoolExecutor

import numpy as np, tifffile
from PIL import Image

from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch

import logging
import argparse
import time
import pickle

from taskqueue import LocalTaskQueue
import igneous.task_creation as tc

def make_info_file(volume_size,resolution,layer_dir,commit=True):
    """ 
    ---PURPOSE---
    Make the cloudvolume info file.
    ---INPUT---
    volume_size     [Nx,Ny,Nz] in voxels, e.g. [2160,2560,1271]
    pix_scale_nm    [size of x pix in nm,size of y pix in nm,size of z pix in nm], e.g. [5000,5000,10000]
    commit          if True, will write the info/provenance file to disk. 
                    if False, just creates it in memory
    """
    info = CloudVolume.create_new_info(
        num_channels = 1,
        layer_type = 'image', # 'image' or 'segmentation'
        data_type = 'uint16', # 
        encoding = 'raw', # other options: 'jpeg', 'compressed_segmentation' (req. uint32 or uint64)
        resolution = resolution, # Size of X,Y,Z pixels in nanometers, 
        voxel_offset = [ 0, 0, 0 ], # values X,Y,Z values in voxels
        chunk_size = [ 1024,1024,1 ], # rechunk of image X,Y,Z in voxels -- only used for downsampling task I think
        volume_size = volume_size, # X,Y,Z size in voxels
        )

    vol = CloudVolume(f'file://{layer_dir}', info=info)
    vol.provenance.description = "Test on spock for profiling precomputed creation"
    vol.provenance.owners = ['zmd@princeton.edu'] # list of contact email addresses
    if commit:
        vol.commit_info() # generates info json file
        vol.commit_provenance() # generates provenance json file
        print("Created CloudVolume info file: ",vol.info_cloudpath)
    return vol

def make_downsample_tasks(vol,mip_start=0,num_mips=3):
    """ 
    ---PURPOSE---
    Make downsamples of the precomputed data
    ---INPUT---
    vol             The cloudvolume.Cloudvolume() object
    mip_start       The mip level to start at with the downsamples
    num_mips        The number of mip levels to create, starting from mip_start
    """
    cloudpath = vol.cloudpath
    tasks = tc.create_downsampling_tasks(
        cloudpath, 
        mip=mip_start, # Start downsampling from this mip level (writes to next level up)
        fill_missing=False, # Ignore missing chunks and fill them with black
        axis='z', 
        num_mips=num_mips, # number of downsamples to produce. Downloaded shape is chunk_size * 2^num_mip
        chunk_size=[ 128, 128, 64 ], # manually set chunk size of next scales, overrides preserve_chunk_size
        preserve_chunk_size=True, # use existing chunk size, don't halve to get more downsamples
      )
    return tasks

def process_slice(z):
    if os.path.exists(os.path.join(progress_dir, str(z))):
        print(f"Slice {z} already processed, skipping ")
        return
    if z > (len(sorted_files) - 1):
        print("Index {z} is larger than (number of slices - 1), skipping")
        return
    print('Processing slice z=',z)
    img_name = sorted_files[z]
    image = Image.open(img_name)
    width, height = image.size 
    array = np.array(image, dtype=np.uint16, order='F')
    array = array.reshape((1, height, width)).T
    vol[:,:, z] = array
    image.close()
    touch(os.path.join(progress_dir, str(z)))
    print("success")
    return

if __name__ == "__main__":
    """ First command line arguments """
    step = sys.argv[1]
    viz_dir = sys.argv[2]
    animal_id = sys.argv[3]
    print(f"Viz_dir: {viz_dir}")
    print(f"Animal id: {animal_id}")
    
    full_sizedatafld = os.path.join('/jukebox/wang/pisano/tracing_output/antero_4x',f'{animal_id}','full_sizedatafld')
    rawdata_path = glob.glob(full_sizedatafld + f'/{animal_id}*647*')[0]
    #zmd added - to grab dims automatically
    img = tifffile.imread(os.path.join(rawdata_path, os.listdir(rawdata_path)[0]))
    layer_name = f'rawdata_{animal_id}'
    layer_dir = os.path.join(viz_dir,layer_name)
    """ Make progress dir """
    progress_dir = mkdir(viz_dir + f'/progress_{layer_name}') # unlike os.mkdir doesn't crash on prexisting 
    """ Figure out volume size in pixels and in nanometers """
    x_dim = img.shape[1]
    y_dim = img.shape[0]
    all_slices = glob.glob(f"{rawdata_path}/*tif") 
    sorted_files = sorted(all_slices)
    z_dim = len(sorted_files)
    
    x_scale_nm, y_scale_nm,z_scale_nm = 10000,10000,10000 

    """ Handle the different steps """
    if step == 'step0':
        print("step 0")
        volume_size = (x_dim,y_dim,z_dim)
        resolution = (x_scale_nm,y_scale_nm,z_scale_nm)
        vol = make_info_file(volume_size=volume_size,layer_dir=layer_dir,resolution=resolution)
    elif step == 'step1':
        print("step 1")
        # Find the individual z planes in the full_sizedatafld - these are the blended images at the raw resolution
        vol = CloudVolume(f'file://{layer_dir}')
        done_files = set([ int(z) for z in os.listdir(progress_dir) ])
        all_files = set(range(vol.bounds.minpt.z, vol.bounds.maxpt.z + 1))

        to_upload = [ int(z) for z in list(all_files.difference(done_files)) ]
        to_upload.sort()
        print(f"Have {len(to_upload)} planes to upload")
        with ProcessPoolExecutor(max_workers=16) as executor:
            executor.map(process_slice, to_upload)

        vol.cache.flush()
    elif step == 'step2': # downsampling
        print("step 2")
        vol = CloudVolume(f'file://{layer_dir}')
        tasks = make_downsample_tasks(vol,mip_start=0,num_mips=4)
        with LocalTaskQueue(parallel=18) as tq:
            tq.insert_all(tasks)


