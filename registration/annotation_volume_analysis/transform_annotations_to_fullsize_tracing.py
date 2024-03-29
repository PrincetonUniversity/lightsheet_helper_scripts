#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 16:15:28 2019

@author: wanglab
"""

import os, numpy as np, sys, time, tifffile, SimpleITK as sitk, multiprocessing as mp
from scipy.ndimage.interpolation import zoom
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.utils.io import makedir, load_memmap_arr, listall, load_kwargs
from tools.registration.register import change_interpolation_order, transformix_command_line_call
from tools.registration.transform_list_of_points import modify_transform_files

def fast_scandir(dirname):
    """ gets all folders recursively """
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders

def resize_helper(iii,zplane,ap0,ap1,ml0,ml1,brain,save_dst):
    """ to upsize annotation volume to raw space """
    arr = zoom(zplane, (ap0/float(ap1), ml0/float(ml1)), order=0)
    tifffile.imsave(os.path.join(save_dst, 
    "%s_annotation_Z%04d.tif" % (os.path.basename(brain),iii)),arr.astype("float32"),compress=6)
    if iii%100==0: print("made z plane # {}".format(iii))
    return
    
if __name__ == "__main__":
    
    #setting paths
    ann = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif" #atlas you registered to!
    scratch_dir = "/jukebox/scratch/zmd/"
    src = "/jukebox/LightSheetTransfer/tp"
    brains = ["PRV_50hr-019", "20201001_10_57_49_hsv_36h_6","20201001_10_01_03_hsv_36h_5",
              "20201001_15_39_26_hsv_28h_4","20201001_17_13_35_hsv_28h_2",
              "20200930_18_34_47_hsv_28hr_3"]
    #cores for parallelization, make sure this matches bash script
    cores = 12
    #channel to use as reference image
    ch = "Ex_642_Em_2"
    #for array job parallelization
    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    # jobid=0
    #set brain name
    brain = os.path.join(src, brains[jobid])
    start = time.time()
    #need to change this config depending on processing pipeline
    cellvol = fast_scandir(os.path.join(brain,ch,"stitched"))[-1] 
    #link to parameters
    a2r0 = os.path.join(brain, "elastix_inverse_transform/TransformParameters.0.txt")
    a2r1 = os.path.join(brain, "elastix_inverse_transform/TransformParameters.1.txt")
    #set destination directory
    braindst = os.path.join(scratch_dir, os.path.basename(brain))
    makedir(braindst)
    aldst = os.path.join(braindst, "transformed_annotations"); makedir(aldst)
    #transformix
    # transformfiles = modify_transform_files(transformfiles=[a2r0, a2r1, r2s0, r2s1], dst = aldst)
    transformfiles = modify_transform_files(transformfiles=[a2r0, a2r1], dst = aldst)
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
    pl0 = sitk.GetArrayFromImage(sitk.ReadImage((os.path.join(cellvol, os.listdir(cellvol)[0]))))
    #get dimensions for proportions
    dv0,ap0,ml0 = len(os.listdir(cellvol)), pl0.shape[0], pl0.shape[1]   
    ml1,ap1,dv1 = tann.shape
    #scale in dv only first and rotate to hor orientation
    bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)
    save_dst = os.path.join(aldst, "single_tifs"); makedir(save_dst)
    #now rotate and scale each in ap and ml
    #parallelize
    iterlst = [(iii,zplane,ap0,ap1,ml0,ml1,brain,save_dst) for iii,zplane in enumerate(bigdvann)]
    p = mp.Pool(cores)
    p.starmap(resize_helper, iterlst)
    p.terminate()
    print("\n\ntook {} minutes for {}\n".format(round((time.time()-start)/60, 2), brain))
