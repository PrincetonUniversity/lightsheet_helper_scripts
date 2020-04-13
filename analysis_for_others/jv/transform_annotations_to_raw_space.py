#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:40:30 2020

@author: wanglab
"""

import os, numpy as np, sys, time
from skimage.external import tifffile
from scipy.ndimage.interpolation import zoom
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.utils.io import makedir
from tools.registration.register import change_interpolation_order, transformix_command_line_call
from tools.registration.transform_list_of_points import modify_transform_files

#setting paths
ann = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
scratch_dir = "/jukebox/scratch/zmd"
src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/processed"

print(sys.argv)
print(os.environ["SLURM_ARRAY_TASK_ID"])
jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

#list of brains
brains = [os.path.join(src, xx) for xx in os.listdir(src)]

#set brain name
brain = brains[jobid]
print("\n**********"+os.path.basename(brain)+"**********\n")

start = time.time()

a2r0 = os.path.join(brain, "clearmap_cluster_output/elastix_cfos_to_auto/TransformParameters.0.txt")
a2r1 = os.path.join(brain, "clearmap_cluster_output/elastix_cfos_to_auto/TransformParameters.1.txt")
r2s0 = os.path.join(brain, "clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters.0.txt")
r2s1 = os.path.join(brain, "clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters.1.txt")

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
    filedata = filedata.replace('(ResultImagePixelType "short")', 
                                '(ResultImagePixelType "float")')
    # Write the file out again
    with open(fl, "w") as file:
      file.write(filedata)
#run transformix
transformix_command_line_call(ann, aldst, transformfiles[-1])

#now zoom out - this is heavy!
transformed_ann = os.path.join(aldst, "result.tif")
tann = tifffile.imread(transformed_ann)
fld = os.path.join(os.path.join(brain, "full_sizedatafld"),
                   os.listdir(os.path.join(brain, "full_sizedatafld"))[0])
pl0 = tifffile.imread(os.path.join(fld, os.listdir(fld)[0]))
dv0,ap0,ml0 = len(os.listdir(fld)), pl0.shape[0], pl0.shape[1]

ml1,ap1,dv1 = tann.shape

#scale in dv only first and rotate to hor orientation
bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)
save_dst = os.path.join(aldst, "single_tifs"); makedir(save_dst)

#now rotate and scale each in ap and ml
for iii,zplane in enumerate(bigdvann):
    arr = zoom(zplane, (ap0/float(ap1), ml0/float(ml1)), order=0)
    tifffile.imsave(os.path.join(save_dst, 
             "%s_annotation_Z%04d.tif" % (os.path.basename(brain), iii)), arr, compress = 6)
    print("\nmade z plane # {} for annotation".format(iii))

#make another folder of only brain volume in raw space
single_tifs = os.listdir(save_dst); single_tifs.sort() #important to keep z-planes consistent
for iii, tif in enumerate(single_tifs):
    img = tifffile.imread(os.path.join(save_dst, tif))
    img[img > 0] = 255
    tifffile.imsave(os.path.join(save_dst,
        "%s_brain_vol_Z%04d.tif" % (os.path.basename(brain), iii)), img.astype("uint8"))
    print("\nmade z plane # {} for brain volume".format(iii))

print("\n\ntook {} minutes for {}\n".format(round((time.time()-start)/60, 2), brain))
