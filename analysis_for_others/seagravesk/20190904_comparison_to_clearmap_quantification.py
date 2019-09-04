#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 11:13:06 2019

@author: wanglab
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from tools.conv_net.utils.io import pairwise_distance_metrics, read_roi_zip

pth = "/jukebox/wang/zahra/kelly_cell_detection_analysis/comparison_to_clearmap/annotated_volumes"

flds = [os.path.join(pth, "sweep/"+xx) for xx in os.listdir(os.path.join(pth, "sweep")) if "max" in xx]
roi_pths = [os.path.join(pth, xx) for xx in os.listdir(pth) if "RoiSet.zip" in xx]; roi_pths.sort() #the sorts here are important, for order of volumes

#these will be zyx 
#note that importing it this way, the z dimension does not start from 0, but neither does clearmaps, so probably ok???
annotated_cells = np.array([np.array([[int(yy) for yy in xx[0].replace(".roi", "").split("-")]for xx in
                            read_roi_zip(roi_pth, include_roi_name=True)]) for roi_pth in roi_pths])

#voxel pair cutoff
cutoffs = [5,7,10]

for cutoff in cutoffs:
    #init dataframe
    print("cutoff: %s \n\n" % cutoff)
    df = pd.DataFrame()
    df["parameters"] = [os.path.basename(xx) for xx in flds]
    df["true_positives"] = np.zeros(len(df))
    df["false_positives"] = np.zeros(len(df))
    df["false_negatives"] = np.zeros(len(df))
    df["f1"] = np.zeros(len(df))
    
    for n, fld in enumerate(flds):
        if n%100 == 0: print(n)
        arr_pths = [os.path.join(fld, xx) for xx in os.listdir(fld) if not os.path.isdir(os.path.join(fld, xx))]; arr_pths.sort() #the sorts here are important, for order of volumes
        detected_cells = np.asarray([np.load(arr_pth) for arr_pth in arr_pths])
        measures = [pairwise_distance_metrics(annotated_cells[i], detected_cells[i], cutoff = cutoff, verbose = False) for i
                    in range(len(arr_pths))] 
        #measures consists of paired, tp, fp, fn, in that order
        tp = sum([measures[i][1] for i in range(len(measures))])
        fp = sum([measures[i][2] for i in range(len(measures))])
        fn = sum([measures[i][3] for i in range(len(measures))])
        precision = tp/(tp+fp); recall = tp/(tp+fn) #calculating precision and recall
        f1 = 2*( (precision*recall)/(precision+recall) ) #calculating f1 score
        df.loc[df.parameters == os.path.basename(fld), "f1"] = f1
        df.loc[df.parameters == os.path.basename(fld), "true_positives"] = tp
        df.loc[df.parameters == os.path.basename(fld), "false_positives"] = fp
        df.loc[df.parameters == os.path.basename(fld), "false_negatives"] = fn
        
    df.to_csv(os.path.join(pth, "clearmap_performance_measures_kelly_cfos_voxelcutoff_%02d_v2.csv" % cutoff))
    


