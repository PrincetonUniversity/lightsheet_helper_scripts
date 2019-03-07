#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 12:40:14 2019

@author: wanglab
"""

import numpy as np, os
from skimage.external import tifffile as tif
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.registration.transform_list_of_points import modify_transform_files
from tools.utils.io import makedir
from tools.registration.transform_list_of_points import point_transformix, unpack_pnts, create_text_file_for_elastix
from tools.registration.register import change_transform_parameter_initial_transform

def generate_transformed_cellcount(points, dst, transformfiles):
    """ makes an input to transformix """
    #set up locations
    transformed_dst = os.path.join(dst, "transformed_points"); makedir(transformed_dst)
    
    #make zyx numpy arry
    zyx = np.asarray(points)
    
    transformfiles = modify_transform_files(transformfiles, transformed_dst) 
    change_transform_parameter_initial_transform(transformfiles[0], "NoInitialTransform")
   
    pretransform_text_file = create_text_file_for_elastix(zyx, transformed_dst)
    
    #run transformix on points
    points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
    
    #convert registered points into structure counts
    converted_points = unpack_pnts(points_file, transformed_dst)   
    
    return converted_points

ann = tif.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif")

dst = "/home/wanglab/Desktop/test1"

pixel_values = np.unique(ann).astype("uint16")

pixel_value = 136 #id that corresponds to structure

structure = np.asarray(np.where(ann == pixel_value)) #gives a tuple of the z,y,x coordinates corresponding to the structure

z = structure[0]; y = structure[1]; x = structure[2]

points = [[z[ii], y[ii], x[ii]] for ii in range(len(z))] #these are the points you would want to map to atlas space

#transform files
r2atl0 = "/home/wanglab/mounts/wang/pisano/tracing_output/antero_4x/20160801_db_cri_02_1200rlow_52hr/elastix/20160801_db_cri_02_1200rlow_52hr_647_017na_1hfds_z7d5um_200msec_10povlp_resized_ch00/sig_to_reg/regtoatlas_TransformParameters.0.txt"
r2atl1 = "/home/wanglab/mounts/wang/pisano/tracing_output/antero_4x/20160801_db_cri_02_1200rlow_52hr/elastix/20160801_db_cri_02_1200rlow_52hr_647_017na_1hfds_z7d5um_200msec_10povlp_resized_ch00/sig_to_reg/regtoatlas_TransformParameters.1.txt"

cell2r0 = "/home/wanglab/mounts/wang/pisano/tracing_output/antero_4x/20160801_db_cri_02_1200rlow_52hr/elastix/20160801_db_cri_02_1200rlow_52hr_647_017na_1hfds_z7d5um_200msec_10povlp_resized_ch00/sig_to_reg/TransformParameters.0.txt"
cell2r1 = "/home/wanglab/mounts/wang/pisano/tracing_output/antero_4x/20160801_db_cri_02_1200rlow_52hr/elastix/20160801_db_cri_02_1200rlow_52hr_647_017na_1hfds_z7d5um_200msec_10povlp_resized_ch00/sig_to_reg/TransformParameters.1.txt"

transformfiles = [r2atl0, r2atl1, cell2r0, cell2r1]
transformfiles = modify_transform_files(transformfiles, dst)

converted_points = generate_transformed_cellcount(points, dst, transformfiles)