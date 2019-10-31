#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:19:28 2019

@author: wanglab
"""

import os, numpy as np, sys, time
from skimage.external import tifffile
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.utils.io import makedir, load_memmap_arr, listall, load_kwargs
from tools.registration.register import change_interpolation_order, transformix_command_line_call
from tools.registration.transform_list_of_points import modify_transform_files
from scipy.ndimage.interpolation import zoom

#setting paths
ann = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif"
scratch_dir = "/jukebox/scratch/zmd"
src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/processed"

#set brain name
brain = os.path.join(src, "an19")

start = time.time()

a2r0 = os.path.join(brain, "clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters.0.txt")
a2r1 = os.path.join(brain, "clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters.1.txt")
r2s0 = os.path.join(brain, "clearmap_cluster_output/elastix_cfos_to_auto/TransformParameters.0.txt")
r2s1 = os.path.join(brain, "clearmap_cluster_output/elastix_cfos_to_auto/TransformParameters.1.txt")

#set destination directory
braindst = os.path.join(scratch_dir, os.path.basename(brain))

makedir(braindst)
    
aldst = os.path.join(braindst, "transformed_annotations"); makedir(aldst)
#
#transformix
transformfiles = modify_transform_files(transformfiles=[a2r0, a2r1, r2s0, r2s1], dst = aldst)
[change_interpolation_order(xx,0) for xx in transformfiles]

#change the parameter in the transform files that outputs 16bit images instead
for fl in transformfiles:# Read in the file
    with open(fl, "r") as file:
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
    # Write the file out again
    with open(fl, "w") as file:
      file.write(filedata)
#run transformix  
transformix_command_line_call(ann, aldst, transformfiles[-1])

#now zoom out - this is heavy!
transformed_ann = os.path.join(aldst, "result.tif")
tann = tifffile.imread(transformed_ann)
flszdt = os.path.join(brain, "full_sizedatafld")
cellfld = [os.path.join(flszdt, xx) for xx in os.listdir(flszdt) if "647" in xx][0]
pl0 = tifffile.imread(os.path.join(cellfld, os.listdir(cellfld)[0]))
dv0,ap0,ml0 = len(os.listdir(cellfld)), pl0.shape[0], pl0.shape[1]

ml1,ap1,dv1 = tann.shape

#scale in dv only first and rotate to hor orientation
bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)
save_dst = os.path.join(aldst, "single_tifs"); makedir(save_dst)

#now rotate and scale each in ap and ml
for iii,zplane in enumerate(bigdvann):
    arr = zoom(zplane, (ap0/float(ap1), ml0/float(ml1)), order=0)
    tifffile.imsave(os.path.join(save_dst, "%s_annotation_Z%04d.tif" % (os.path.basename(brain), iii)), arr, compress = 6)
    print("\nmade z plane # {}".format(iii))
          
print("\n\ntook {} minutes for {}\n".format(round((time.time()-start)/60, 2), brain))
    