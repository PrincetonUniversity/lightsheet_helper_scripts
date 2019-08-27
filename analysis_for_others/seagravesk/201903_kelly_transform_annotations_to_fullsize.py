#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 16:41:30 2019

@author: wanglab
"""

import os, numpy as np, sys, time
from skimage.external import tifffile
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.utils.io import makedir, load_memmap_arr, listall, load_kwargs
from tools.registration.register import change_interpolation_order, transformix_command_line_call
from tools.registration.transform_list_of_points import modify_transform_files
from tools.registration.transform_cell_counts import get_fullsizedims_from_kwargs
from scipy.ndimage.interpolation import zoom

#setting paths
ann = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans.nrrd"
scratch_dir = "/jukebox/scratch/kellyms"
src = "/jukebox/wang/seagravesk/lightsheet/201908_cfos_left_side_only_registration"

print(sys.argv)
print(os.environ["SLURM_ARRAY_TASK_ID"])
jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

#list of brains    
brains = [os.path.join(src, xx) for xx in os.listdir(src)]

#set brain name
brain = brains[jobid]

start = time.time()

kwargs = load_kwargs(brain)

#accessing parameter dictionary
cellvol = [xx for xx in kwargs["volumes"] if xx.ch_type == "cellch"][0]

a2r0 = [xx for xx in listall(cellvol.inverse_elastixfld) if "atlas2reg_TransformParameters.0" in xx and "cellch" in xx][0]
a2r1 = [xx for xx in listall(cellvol.inverse_elastixfld) if "atlas2reg_TransformParameters.1" in xx and "cellch" in xx][0]
r2s0 = [xx for xx in listall(cellvol.inverse_elastixfld) if "reg2sig_TransformParameters.0" in xx and "cellch" in xx][0]
r2s1 = [xx for xx in listall(cellvol.inverse_elastixfld) if "reg2sig_TransformParameters.1" in xx and "cellch" in xx][0]

#set destination directory
braindst = os.path.join(scratch_dir, os.path.basename(brain)+"_dorsal")
makedir(braindst)
    
aldst = os.path.join(braindst, "transformed_annotations"); makedir(aldst)

#transformix
transformfiles = modify_transform_files(transformfiles=[a2r0, a2r1, r2s0, r2s1], dst = aldst)
[change_interpolation_order(xx,0) for xx in transformfiles]
transformix_command_line_call(ann, aldst, transformfiles[-1])

#now zoom out - this is heavy!
transformed_ann = os.path.join(aldst, "result.tif")
tann = tifffile.imread(transformed_ann)
dv0,ap0,ml0 = get_fullsizedims_from_kwargs(kwargs)

ml1,ap1,dv1 = tann.shape

#scale in dv only first and rotate to hor orientation
bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)

#make memmap
annout = os.path.join(aldst, "transformed_annotations.npy")
annout = load_memmap_arr(annout, mode="w+", dtype="float32", shape=(dv0,ap0,ml0))

#now rotate and scale each in ap and ml
for iii,zplane in enumerate(bigdvann):
    annout[iii] = zoom(zplane, (ap0/float(ap1), ml0/float(ml1)), order=0)
    if iii%50==0: 
        annout.flush()
        print(iii)
annout.flush()

print("\n\ntook {} minutes to make merged maps for {}\n".format(round((time.time()-start)/60, 2), brain))
    