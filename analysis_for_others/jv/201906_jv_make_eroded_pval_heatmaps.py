#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:49:32 2019

@author: wanglab
"""

import numpy as np, os, tifffile
import multiprocessing as mp
from scipy.ndimage.morphology import distance_transform_edt
from skimage.filters import gaussian
os.chdir("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.Analysis.Voxelization import voxelize
import math

def make_heatmaps(pth):
    """ 
    makes clearmap style heatmaps for eroded data
    """
    
    #print brainname being analyzed
    print(os.path.basename(pth))
    
    #set destination for heatmap
    dst = os.path.join(pth, "cells_heatmap_60um_erosion.tif")
    
    #only make if not available
#    if not os.path.exists(dst):
    #set output folder
    output = os.path.join(pth, "clearmap_cluster_output")
    
    points = np.load(os.path.join(output, "cells_transformed_to_Atlas.npy"))
    
    #make an empty array to put cells on
    cell_map = np.ones((352, 640, 540)) #dims of horizontal atlas
    
    #make a cell map
    for z,y,x in points:
        try:
            cell_map[int(math.ceil(z)), int(math.ceil(y)), int(math.ceil(x))] = 45000
        except Exception as e:
            print(e)
            
    #erode the cell map
    #erode maps like you would do annotation file
    #eroding edges AND ventricles
    zyx_scale = (20, 20, 20)
    microns_to_erode = 60
            
    #NOTE THIS ESSENTIALLY SCALES PIXEL SPACE*****
    distance_space_inside = distance_transform_edt(cell_map.astype('bool'), sampling = zyx_scale)*-1 #INSIDE
    distance_space_inside = np.abs(distance_space_inside)
    mask = np.copy(distance_space_inside)
    mask[distance_space_inside <= microns_to_erode] = 0
    
    #zero out edges
    cell_map_e = np.copy(cell_map)
    cell_map_e[mask==0]=0
    
    #now get the eroded points that were transformed on the cell map
    z,y,x = np.where(cell_map_e > 1)
    
    #remake points into 3d array
    points_e = np.asarray([[x0,y0,z0] for z0, y0, x0 in zip(z,y,x)])
    
    #run clearmap style blurring - this will be horizontal
    vox_params = {"method": "Spherical", "size": (15, 15, 15), "weights": None}
    #transform to SAGITTAL
    vox = voxelize(np.asarray(points_e), dataSize = (540, 640, 352), **vox_params)

    tifffile.imsave(dst, vox.astype("int32"))
            
    return dst

if __name__ == "__main__":
    
    src = "/home/wanglab/mounts/wang/Jess/lightsheet_output/201904_ymaze_cfos/processed"
    
    brains = [os.path.join(src, xx) for xx in os.listdir(src)]
    
    p = mp.Pool(12)
    p.map(make_heatmaps, brains)