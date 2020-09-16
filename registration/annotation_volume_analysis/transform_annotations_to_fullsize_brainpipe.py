#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 18:03:50 2020

@author: wanglab
"""

import os, numpy as np, sys, time
import tifffile, SimpleITK as sitk
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.utils.io import makedir
from tools.registration.register import change_interpolation_order, transformix_command_line_call
from tools.registration.transform_list_of_points import modify_transform_files
from scipy.ndimage.interpolation import zoom

if __name__ == "__main__":
    
    #setting paths
    ann = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
    scratch_dir = "/jukebox/scratch/zmd/"
    src = "/jukebox/wang/pisano/tracing_output/bl6_ts"
    brains = ["20150804_tp_bl6_ts04"]
    #for array job parallelization
    # print(os.environ["SLURM_ARRAY_TASK_ID"])
    # jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid=0
    #set brain name
    brain = os.path.join(src, brains[jobid])

    start = time.time()
    
    #accessing parameter dictionary
    cellvol = "/jukebox/wang/pisano/tracing_output/bl6_ts/20150804_tp_bl6_ts04/full_sizedatafld/20150804_tp_bl6_ts04_488w_647_z3um_250msec_1hfds_ch01"
    
    a2r0 = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125/20150804_tp_bl6_ts04/elastix_inverse_transform/TransformParameters.0.txt"
    a2r1 = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125/20150804_tp_bl6_ts04/elastix_inverse_transform/TransformParameters.1.txt"
    # r2s0 = "/jukebox/wang/pisano/tracing_output/bl6_ts/20150804_tp_bl6_ts04/elastix_inverse_transform/cellch_20150804_tp_bl6_ts04_555_z3um_70msec_3hfds/20150804_tp_bl6_ts04_555_z3um_70msec_3hfds_resized_ch00_resampledforelastix_atlas2reg2sig/reg2sig_TransformParameters.0.txt"
    # r2s1 = "/jukebox/LightSheetTransfer/tp/20200701_12_55_28_20170207_db_bl6_crii_rpv_01/elastix_inverse_transform/TransformParameters.1.txt"
    
    #set destination directory
    braindst = os.path.join(scratch_dir, os.path.basename(brain))
    
    makedir(braindst)
        
    aldst = os.path.join(braindst, "transformed_annotations"); makedir(aldst)
    
    #transformix
    # transformfiles = modify_transform_files(transformfiles=[a2r0, a2r1, r2s0], dst = aldst)
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
    
    # #now zoom out - this is heavy!
    transformed_ann = os.path.join(aldst, "result.tif")
    tann = tifffile.imread(transformed_ann).astype("float32")
    pl0 = sitk.GetArrayFromImage(sitk.ReadImage((os.path.join(cellvol, os.listdir(cellvol)[0]))))
    dv0,ap0,ml0 = len(os.listdir(cellvol)), pl0.shape[0], pl0.shape[1]
    
    ml1,ap1,dv1 = tann.shape
    
    #scale in dv only first and rotate to hor orientation
    bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)
    save_dst = os.path.join(aldst, "single_tifs"); makedir(save_dst)
    
    #now rotate and scale each in ap and ml
    for iii,zplane in enumerate(bigdvann):
        arr = zoom(zplane, (ap0/float(ap1), ml0/float(ml1)), order=0)
        tifffile.imsave(os.path.join(save_dst, "%s_annotation_Z%04d.tif" % (os.path.basename(brain), iii)), arr.astype("float32"), compress = 6)
        if iii%100==0: print("\nmade z plane # {}".format(iii))
              
    print("\n\ntook {} minutes for {}\n".format(round((time.time()-start)/60, 2), brain))