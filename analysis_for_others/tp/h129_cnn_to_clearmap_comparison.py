#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:40:40 2019

@author: wanglab
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from tools.conv_net.utils.io import pairwise_distance_metrics, read_roi_zip

pth = "/home/wanglab/Desktop/cnn_to_clearmap/all_volumes"
dst = os.path.join(pth, "results")
if not os.path.exists(dst): os.mkdir(dst)

src = os.path.join(pth, "cell_arrays")
brains = [os.path.join(src, xx) for xx in os.listdir(src)] #only using the ones from the new sweep
src = "/jukebox/wang/pisano/conv_net/annotations/all_better_res/h129/input_files/"
vols = [os.path.join(src, os.path.basename(fl)+".tif") for fl in brains]
roipths = [os.path.join(src, os.path.basename(fl)+".roi.zip") for fl in brains]
#these will be zyx 
#note that importing it this way, the z dimension does not start from 0, but neither does clearmaps, so probably ok???
annotated_cells = np.array([np.array([[int(yy) for yy in xx[0].replace(".roi", "").split("-")]for xx in
                            read_roi_zip(roipth, include_roi_name=True)]) for roipth in roipths])

#voxel pair cutoff
cutoff = 30

for br, brain in enumerate(brains):
    fls = [os.path.join(brain, xx) for xx in os.listdir(brain)]; fls.sort()
    #init dataframe
    print("cutoff: %s \n\n" % cutoff)
    df = pd.DataFrame()
    df["parameters"] = [os.path.basename(xx) for xx in fls]
    
    df["tp"] = np.zeros(len(df))
    df["fp"] = np.zeros(len(df))
    df["fn"] = np.zeros(len(df))
    df["f1"] = np.zeros(len(df))
    
    for n, fl in enumerate(fls):
        if n%100 == 0: print(n)
        detected_cells = np.load(fl)
        paired, tp, fp, fn = pairwise_distance_metrics(annotated_cells[br], detected_cells, cutoff = cutoff, verbose = False) 
        
        try: 
            precision = tp/(tp+fp)
            recall = tp/(tp+fn) #calculating precision and recall
            f1 = 2*( (precision*recall)/(precision+recall) ) #calculating f1 score
        except Exception as e:
            print(e)
            f1 = np.nan #if tp, fn, etc. are 0                
            
        df.loc[df.parameters == os.path.basename(fl), "f1"] = f1
        df.loc[df.parameters == os.path.basename(fl), "tp"] = tp
        df.loc[df.parameters == os.path.basename(fl), "fp"] = fp
        df.loc[df.parameters == os.path.basename(fl), "fn"] = fn

    #export csv per brain/volume                
    df.to_csv(os.path.join(dst, "%s.csv" % os.path.basename(brain)))
