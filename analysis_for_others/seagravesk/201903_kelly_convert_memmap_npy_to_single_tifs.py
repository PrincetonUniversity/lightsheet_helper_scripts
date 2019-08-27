#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:52:03 2019

@author: wanglab
"""

import numpy as np, os, sys
from skimage.external import tifffile

print(sys.argv)
print(os.environ["SLURM_ARRAY_TASK_ID"])
jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

src = "/jukebox/scratch/kellyms"

pths = [os.path.join(src, xx) for xx in os.listdir(src) if "dorsal" in xx or "ventral" in xx]

if jobid > len(pths)+1:
    print("array jobs greater than number of brains")
else:
    print("processing brain: {}".format(pths[jobid]))
    
    arrpth = os.path.join(pths[jobid], "transformed_annotations/transformed_annotations.npy")
    dst = os.path.join(pths[jobid], "annotations_as_single_tifs")
    if not os.path.exists(dst): 
        os.mkdir(dst)
        arr = np.lib.format.open_memmap(arrpth, dtype = "float32", mode = "r")
        
        print("\nfile with shape: {}".format(arr.shape))
        
        for z in range(arr.shape[0]):
            tifffile.imsave(os.path.join(dst, "annotation_Z{}.tif".format(str(z).zfill(4))), arr[z], compress = 6)
            print("\nmade z plane # {}".format(z))