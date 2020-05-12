#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 14:42:33 2020

@author: wanglab
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from tools.conv_net.utils.io import pairwise_distance_metrics, read_roi_zip

if __name__ == "__main__":

    #paths and imports
    src = "/jukebox/wang/zahra/kelly_cell_detection_analysis/comparison_to_clearmap/may2020/"
    dst = os.path.join(src, "results")
    brains = [os.path.join(src, xx) for xx in os.listdir(src) if "790" in xx]; brains.sort()
    
    imgsrc = "/jukebox/wang/zahra/kelly_cell_detection_analysis/comparison_to_clearmap/annotated_volumes"
    vols = [os.path.join(imgsrc, os.path.basename(fl)+".tif") for fl in brains]; vols.sort()
    roipths = [os.path.join(imgsrc, xx) for xx in os.listdir(imgsrc) if "RoiSet.zip" in xx]; roipths.sort()
    #these will be zyx 
    #note that importing it this way, the z dimension does not start from 0, but neither does clearmaps, so probably ok???
    annotated_cells = np.array([np.array([[int(yy) for yy in xx.replace(".roi", "").split("-")] for xx in
                                read_roi_zip(roipth, points=True)]) for roipth in roipths])
    
    #voxel pair cutoff
    cutoff = 10
    
    #make dest
    if not os.path.exists(dst): os.mkdir(dst)
    
    for br, brain in enumerate(brains):
        fls = [os.path.join(brain, xx) for xx in os.listdir(brain)]; fls.sort()
        #init dataframe
        print("brain: %s\ncutoff: %s \n\n" % (os.path.basename(brain), cutoff))
        df = pd.DataFrame()
        df["parameters"] = [os.path.basename(xx) for xx in fls]
        
        df["tp"] = np.zeros(len(df))
        df["fp"] = np.zeros(len(df))
        df["fn"] = np.zeros(len(df))
        df["f1"] = np.zeros(len(df))
        
        for n, fl in enumerate(fls):
            if n%50 == 0: print(n)
            detected_cells = np.load(fl, allow_pickle=True)
            paired, tp, fp, fn = pairwise_distance_metrics(annotated_cells[br], 
                        detected_cells, cutoff = cutoff, verbose = False) 
            
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