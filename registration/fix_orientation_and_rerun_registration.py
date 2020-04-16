#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 17:28:49 2020

@author: wanglab
"""

import tifffile, numpy as np, matplotlib.pyplot as plt, os

pth = "/jukebox/LightSheetData/falkner-mouse/scooter/clearmap_processed"

brains = [os.path.join(pth, xx) for xx in os.listdir(pth)]

for brain in brains:
    print("\n"+os.path.basename(brain)+"\n")
    tifs = [os.path.join(brain+"/clearmap_cluster_output", xx) for xx 
            in os.listdir(brain+"/clearmap_cluster_output") if "tif" in xx]
    for tif in tifs:
        print("\n*********************"+os.path.basename(tif)+"*********************\n")
        arr = tifffile.imread(tif)
        arrflip = np.flip(arr, 0)
        tifffile.imsave(tif, arrflip.astype("uint16")) #overwrite file
        print("\nREORIENTATION SUCCESFULL!\n")