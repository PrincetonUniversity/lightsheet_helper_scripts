#! /bin/env python

import os, sys
import glob
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from PIL import Image

from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch

import logging
import argparse
import time
import pickle

from taskqueue import LocalTaskQueue
import igneous.task_creation as tc
from precomputed_utils import calculate_chunks, calculate_factors

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
    vol.provenance.owners = ['ahoag@princeton.edu'] # list of contact email addresses
    if commit:
        vol.commit_info() # generates info json file
        vol.commit_provenance() # generates provenance json file
        print("Created CloudVolume info file: ",vol.info_cloudpath)
    return vol

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
    return "success"

if __name__ == "__main__":
    """ First command line arguments """
    step = sys.argv[1]
    viz_dir = sys.argv[2]
    image_resolution = sys.argv[3]
    channel_name = sys.argv[4]
    cpus = os.cpu_count()

    """ Load param dictionary """
    param_file = viz_dir + f'/precomputed_params_{image_resolution}_ch{channel_name}.p'
    with open(param_file,'rb') as pkl_file:
        param_dict = pickle.load(pkl_file)
    blended_data_path = param_dict['blended_data_path']
    print("Looking for blended data in:")
    print(blended_data_path)

    for dirpath, dirnames, filenames in os.walk(blended_data_path):
        if not dirnames:
            slices_dir = dirpath
            break
    layer_name = param_dict['layer_name']
    if image_resolution != "3.6x" and image_resolution != "4x":
        sys.exit(f"Image resolution must be 3.6x or 4x. Instead it is {image_resolution} which is not supported")
    z_scale_nm = int(float(param_dict['z_step'])*1000) # to convert from microns to nm. Most likely 2 microns for 3.6x
    
    # Make directories for orig layer, destination layer 
    # orig - just for uploading mip=-1
    orig_layer_name = layer_name + '_rechunkme'
    orig_layer_dir = os.path.join(viz_dir,orig_layer_name)
    mkdir(orig_layer_dir)
    progress_dir = mkdir(viz_dir + f'/progress_{orig_layer_name}') # unlike os.mkdir doesn't crash on prexisting 

    # dest - where the rechunked layer will live
    dest_layer_dir = os.path.join(viz_dir,layer_name)
    mkdir(dest_layer_dir)
    rechunked_cloudpath = f'file://{dest_layer_dir}'
    # Figure out volume size in pixels and in nanometers
    all_slices = glob.glob(f"{slices_dir}/*tif")  
    assert len(all_slices) > 0
    random_slice = all_slices[0]
    random_im = Image.open(random_slice)
    x_dim,y_dim = random_im.size
    random_im.close()
    z_dim = len(all_slices)    
    if image_resolution == "3.6x":
        x_scale_nm, y_scale_nm = 1866,1866 
    elif image_resolution == "4x":
        x_scale_nm, y_scale_nm = 1630,1630

    # Handle the different steps 
    if step == 'step0':
        print("step 0, making info file")
        volume_size = (x_dim,y_dim,z_dim)
        resolution = (x_scale_nm,y_scale_nm,z_scale_nm)
        vol = make_info_file(volume_size=volume_size,
            layer_dir=orig_layer_dir,
            resolution=resolution)
    elif step == 'step1':
        print("step 1, making full resolution layer at orig chunk size")
        sorted_files = sorted(all_slices)
        vol = CloudVolume(f'file://{orig_layer_dir}')
        done_files = set([ int(z) for z in os.listdir(progress_dir) ])
        all_files = set(range(vol.bounds.minpt.z, vol.bounds.maxpt.z))

        to_upload = [ int(z) for z in list(all_files.difference(done_files)) ]
        to_upload.sort()
        print(f"Have {len(to_upload)} planes to upload")
        with ProcessPoolExecutor(max_workers=cpus) as executor:
            for job in executor.map(process_slice,to_upload):
                try:
                    print(job)
                except Exception as exc:
                    print(f'generated an exception: {exc}')

    elif step == 'step2': # transfer tasks
        orig_vol = CloudVolume(f'file://{orig_layer_dir}')
        # makes a new layer with mip=0 chunks of size: [128,128,64] 
        first_chunk = calculate_chunks(downsample='full',mip=0) 
        tq = LocalTaskQueue(parallel=cpus)

        tasks = tc.create_transfer_tasks(orig_vol.cloudpath, dest_layer_path=rechunked_cloudpath, 
            chunk_size=first_chunk, mip=0, skip_downsamples=True)
        print(len(tasks))
        tq.insert(tasks)
        tq.execute()

    elif step == 'step3': # downsampling
        print("step 3, downsampling")
        tq = LocalTaskQueue(parallel=cpus)
        downsample="full"
        mips = [0,1,2,3,4] # At each mip level, create the next mip level. That's why 0 is in the list
        for mip in mips:
            print(f"Mip: {mip}")
            cv = CloudVolume(rechunked_cloudpath, mip)
            chunks = calculate_chunks(downsample, mip)
            factors = calculate_factors(downsample, mip)
            print(f"Chunk size: {chunks}")
            print(f"Downsample factors: {factors}")
            tasks = tc.create_downsampling_tasks(cv.layer_cloudpath, 
                mip=mip, num_mips=1, factor=factors, preserve_chunk_size=False,
                compress=True, chunk_size=chunks)
            tq.insert(tasks)
            tq.execute()
            print()

