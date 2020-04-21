#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 12:21:39 2020

@author: wanglab
"""

import os, numpy as np, tifffile, sys, matplotlib.pyplot as plt
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.imageprocessing.orientation import fix_orientation
sys.path.append("/jukebox/wang/zahra/python/ClearMapCluster")
from ClearMap.cluster.utils import load_kwargs
from ClearMap.cluster.par_tools import output_analysis
      
if __name__ == "__main__":
    
    src = "/jukebox/LightSheetData/falkner-mouse/scooter/clearmap_processed"
    
    brains = [os.path.join(src, fld) for fld in os.listdir(src)]
    
    for brain in brains:
        print("\n"+os.path.basename(brain)+"\n")
        cells_pth = os.path.join(brain, "clearmap_cluster_output/cells-allpoints.npy")
        arr = np.load(cells_pth)
        fszfld = os.path.join(os.path.join(brain, "full_sizedatafld"),
                  [xx for xx in os.listdir(os.path.join(brain, 
                "full_sizedatafld")) if "647" in xx][0])
        ydim, xdim = tifffile.imread(os.path.join(fszfld, os.listdir(fszfld)[0])).shape
        zdim = len(os.listdir(fszfld))
        cell_map = np.zeros((zdim,ydim,xdim), dtype = "bool") 
        for x,y,z in arr:
            try:
                cell_map[z,y,x] = True
            except Exception as e:
                print(e)
                
        #flip the z axis in cells
        cells_flip = fix_orientation(cell_map, axes = ("0", "1", "-2"))
        arr_flip = np.where(cells_flip == True)
        #convert to x,y,z 3d array
        arr_flip = np.array([np.array([arr_flip[2][i],arr_flip[1][i],
                                        arr_flip[0][i]]) for i in range(len(arr_flip[0]))])
        #save out/overwrite
        np.save(cells_pth, arr_flip)
        #re-run step 6 to threshold and transform points
        params = load_kwargs(brain)
        output_analysis(threshold = (750, 8000), row = (2,2), 
                        check_cell_detection = False, **params) 
        