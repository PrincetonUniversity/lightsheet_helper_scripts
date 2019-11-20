#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:34:18 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, cv2, tifffile as tif
from skimage.morphology import ball

#first, make a map of cells
pth = "/jukebox/wang/pisano/tracing_output/antero_4x/20170115_tp_bl6_lob6b_ml_04/3dunet_output/pooled_cell_measures/20170115_tp_bl6_lob6b_ml_04_cell_measures.csv"
svpth = "/jukebox/scratch/zmd/20170115_tp_bl6_lob6b_ml_04/cell_map_single_tifs"
if not os.path.exists(svpth): os.mkdir(svpth)

df = pd.read_csv(pth)

#make into array
zyx = np.squeeze(np.array([(df["z"].values, df["y"].values, df["x"].values)])).T #cells are counted in horizontal volumes

fszdt = "/jukebox/wang/pisano/tracing_output/antero_4x/20170115_tp_bl6_lob6b_ml_04/full_sizedatafld/20170115_tp_bl6_lob6b_ml_04_4x_647_008na_1hfds_z7d5um_75msec_10povlp_ch00"
pln0 = tif.imread(os.path.join(fszdt, os.listdir(fszdt)[0]))
dims = (pln0.shape[0], pln0.shape[1])

pln = np.unique(zyx[:,0]) #unique planes
#fill map
for z in range(len(os.listdir(fszdt))):
    if z%10 == 0: print(z)
    
    cell_map = np.zeros(dims)
    zyx_per_pln = np.array([xx for xx in zyx if xx[0] == z])
    
    #if cells exists in this plane
    if not len(zyx_per_pln) == 0:
        for z,y,x in zyx_per_pln:
            try:
                cell_map[y,x] = 255 #z dilation
            except Exception as e:
                print(e)
            #apply x y dilation
        r = 8
        selem = ball(r)[int(r/2)]
        # cell_map = cell_map.astype("uint8")
        cell_map = cv2.dilate(cell_map, selem, iterations = 1)
    
    #else just save the 0 plane
    tif.imsave(os.path.join(svpth, "20170115_tp_bl6_lob6b_ml_04_cell_map_Z%04d.tif" % z), cell_map.astype("uint8"))
