#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:40:40 2019

@author: wanglab
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from tools.conv_net.utils.io import pairwise_distance_metrics, read_roi_zip

pth = "/home/wanglab/Desktop/cnn_to_clearmap"
dst = os.path.join(pth, "results")
if not os.path.exists(dst): os.mkdir(dst)

fls = [os.path.join(pth, xx) for xx in os.listdir(pth)] #only using the ones from the new sweep
vol = "/jukebox/wang/pisano/conv_net/annotations/all_better_res/h129/input_files/20170204_tp_bl6_cri_1000r_02_1hfds_647_0010na_25msec_z7d5um_10povlap_ch00_z200-400_y1000-1350_x2050-2400.tif"
roipth = "/jukebox/wang/pisano/conv_net/annotations/all_better_res/h129/input_files/20170204_tp_bl6_cri_1000r_02_1hfds_647_0010na_25msec_z7d5um_10povlap_ch00_z200-400_y1000-1350_x2050-2400.roi.zip"
#these will be zyx 
#note that importing it this way, the z dimension does not start from 0, but neither does clearmaps, so probably ok???
annotated_cells = np.array([[int(yy) for yy in xx[0].replace(".roi", "").split("-")]for xx in
                            read_roi_zip(roipth, include_roi_name=True)])

#voxel pair cutoff
cutoffs = [30]

for cutoff in cutoffs:
    #init dataframe
    print("cutoff: %s \n\n" % cutoff)
    df = pd.DataFrame()
    df["parameters"] = [os.path.basename(xx) for xx in fls]
    
    df["tp"] = np.zeros(len(df))
    df["fp"] = np.zeros(len(df))
    df["fn"] = np.zeros(len(df))
    df["f1"] = np.zeros(len(df))
    
    for n, fld in enumerate(fls):
        if n%100 == 0: print(n)
        detected_cells = np.load(fld)
        paired, tp, fp, fn = pairwise_distance_metrics(annotated_cells, detected_cells, cutoff = cutoff, verbose = False)
        try: 
            precision = tp/(tp+fp)
            recall = tp/(tp+fn) #calculating precision and recall
            f1 = 2*( (precision*recall)/(precision+recall) ) #calculating f1 score
        except Exception as e:
            print(e)
            f1 = np.nan #if tp, fn, etc. are 0    
            
        df.loc[df.parameters == os.path.basename(fld), "f1"] = f1
        df.loc[df.parameters == os.path.basename(fld), "tp"] = tp
        df.loc[df.parameters == os.path.basename(fld), "fp"] = fp
        df.loc[df.parameters == os.path.basename(fld), "fn"] = fn
                
    df.to_csv(os.path.join(dst, "clearmap_performance_measures_tp_cutoff30.csv"))
