#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:55:49 2019

@author: wanglab
"""

import os, numpy as np, time
from skimage.external import tifffile
import matplotlib.gridspec as gridspec
from scipy.ndimage.interpolation import zoom
from collections import Counter
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.utils.io import makedir, load_dictionary, load_np, listdirfull
from tools.registration.transform_list_of_points import create_text_file_for_elastix, point_transformix, modify_transform_files, unpack_pnts
from tools.registration.register import change_transform_parameter_initial_transform
from tools.registration.transform_cell_counts import points_resample, get_fullsizedims_from_kwargs
os.chdir("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.cluster.utils import load_kwargs


def generate_transformed_cellcount(dataframe, dst, transformfiles, dct, verbose=False):
    """Function to take a csv file and generate an input to transformix
    
    Inputs
    ----------------
    dataframe = preloaded pandas dataframe
    dst = destination to save files
    transformfiles = list of all elastix transform files used, and in order of the original transform****
    lightsheet_parameter_file = .p file generated from lightsheet package
    """
    #set up locations
    transformed_dst = os.path.join(dst, "transformed_points"); makedir(transformed_dst)
    
    #make zyx numpy arry
    zyx = dataframe
    
    #adjust for reorientation THEN rescaling, remember full size data needs dimension change releative to resample
    kwargs = load_dictionary(dct)
    zyx = zyx[:, [2, 1, 0]] #fix orientation  
    fullsizedimensions = get_fullsizedims_from_kwargs(kwargs) #don"t get from kwargs["volumes"][0].fullsizedimensions it"s bad! use this instead
    zyx = points_resample(zyx, original_dims = (fullsizedimensions[2], fullsizedimensions[1], fullsizedimensions[0]), 
                          resample_dims = tifffile.imread(os.path.join(fld, "clearmap_cluster_output/cfos_resampled.tif")).shape, 
                          verbose = verbose)[:, :3]
   
    #make into transformix-friendly text file
    pretransform_text_file = create_text_file_for_elastix(zyx, transformed_dst)
        
    #copy over elastix files
    transformfiles = modify_transform_files(transformfiles, transformed_dst) 
    change_transform_parameter_initial_transform(transformfiles[0], "NoInitialTransform")
   
    #run transformix on points
    points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
    
    #convert registered points into structure counts
    converted_points = unpack_pnts(points_file, transformed_dst)   
    
    return converted_points
#%%
if __name__ == "__main__":

    #set up
    dst = "/jukebox/wang/Jess/lightsheet_output/201812_development/forebrain/qc"; makedir(dst)
    #make logs dir
    makedir(os.path.join(dst, "errors"))
    lst = [xx for xx in listdirfull("/jukebox/wang/Jess/lightsheet_output/201812_development/forebrain/processed")]
    transform = "all";#both for regwatlas, and only affine for sig adn reg #"all", "single": don"t consider reg with sig at all
    verbose = True
    atl = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
    #loop
    for fld in lst:
        start = time.time()
        print(fld)
        kwargs = load_kwargs(fld)        
        
        #####check cell transform
        #clearmap transformed to atlas cell dataframe
        dataframe = load_np(os.path.join(fld, "clearmap_cluster_output/cells_transformed_to_Atlas.npy"))
        
        #load and convert to single voxel loc
        zyx = np.asarray([str((int(xx[2]), int(xx[1]), int(xx[0]))) for xx in load_np(dataframe)])
        zyx_cnt = Counter(zyx)
                            
        #atlas
        atl = tifffile.imread(atl)
        atl_cnn = np.zeros_like(atl)
        errors = []
        for zyx,v in zyx_cnt.iteritems():
            z,y,x = [int(xx) for xx in zyx.replace("(","",).replace(")","").split(",")]
            try:
                atl_cnn[z,y,x] = v
            except Exception, e:
                print e
                errors.append(e)
        if len(errors)>0:
            with open(os.path.join(dst, "errors/{}_errors.txt".format(os.path.basename(fld))), "a") as flll:
                for err in errors:
                    flll.write(str(err)+"\n")
                flll.close()
        merged = np.stack([atl, atl_cnn, np.zeros_like(atl)], -1)
        #reorient to horizontal
        merged = np.swapaxes(merged, 0, 2)
        tifffile.imsave(os.path.join(dst, "{}_points_merged.tif".format(os.path.basename(fld))), merged)
        print("\ntook {} seconds to make merged maps\n".format(time.time()-start))