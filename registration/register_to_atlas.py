#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:09:03 2018

@author: wanglab
"""

import sys, numpy as np
from skimage.external import tifffile

#load atlas and annotations file
ann = tifffile.imread('/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif') #load annotations

atlas = tifffile.imread('/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif') #load atlas


def zero_atlas(ann, atlas, straight_cut = False, verbose = True):
    '''Function to zero out piece of atlas that doesn't exist in the scans. Good for cut brains imaged under lightsheet.
    Keep in mind that 'titled cut' only works if the area needed to be zeroed out is a square; otherwise see mask/dilation functionality in 
    tools.registration in lightsheet package
    Inputs:
        ann = annotation file, defaults to pma atlas
        atlas = atlas file
    
    Returns:
        path to ann and atlas file
    '''
    #find atlas and annotations file shape - triple of z,y,x coordinates
    ann_shape = ann.shape    
    atlas_shape = atlas.shape    
    
    #straight cut
    #'zero' out cerebellum that doesn't exist in the scans
    if straight_cut:
        if verbose: sys.stdout.write('Zeroing the atlas below the midbrain...'); sys.stdout.flush()
        ann[:,390:,:] = np.zeros((456, 138, 320))
        atlas[:,390:,:] = np.zeros((456, 138, 320))
    
    #tilted cut - 201810 cfos data
    if not straight_cut:
        #set points for triangle
        y1 = 461; y2 = 411; x1 = 0; x2 = 115
        if verbose: sys.stdout.write('\n    Using the coordinates: \n\
                                     y1: {}, y2: {}, x1: {}, x2: {} \n\
                                     y1 > y2, x2 > x1'.format(y1, y2, x1, x2)); sys.stdout.flush()
        
        #annotations
        ann_arr = ann[:, y2:y1, x1:x2] #z,y,x convention
        ann_arr_tril = np.flip(np.tril(np.flip(ann_arr, axis=0)), axis=0) #to zero out on the opposite of the identity line
        ann[:, y2:y1, x1:x2] = ann_arr_tril
        ann[:, y1:ann_shape[1], :] = np.zeros((ann_shape[0], ann_shape[1]-y1, ann_shape[2]))
        
        #atlas
        atlas_arr = atlas[:, y2:y1, x1:x2] #z,y,x convention
        atlas_arr_tril = np.flip(np.tril(np.flip(atlas_arr, axis=0)), axis=0)
        atlas[:, y2:y1, x1:x2] = atlas_arr_tril
        atlas[:, y1:atlas_shape[1], :] = np.zeros((ann_shape[0], atlas_shape[1]-y1, atlas_shape[2]))
    
        #save out tifs in the appropriate location
        tifffile.imsave('/jukebox/LightSheetTransfer/atlas/201810_jess_cfos/annotation_template_25_sagittal_forDVscans_forebrain_only_modified.tif', ann)
        tifffile.imsave('/jukebox/LightSheetTransfer/atlas/201810_jess_cfos/average_template_25_sagittal_forDVscans_forebrain_only_modified.tif', atlas)