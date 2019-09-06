#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 17:15:41 2019

@author: wanglab
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt, pickle
from tools.conv_net.utils.io import pairwise_distance_metrics, read_roi_zip

src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/clearmap_accuracy_quantification"
dst = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/clearmap_accuracy_quantification/results"
#keeping this order bc that is how the rois are save in the pickle file
rois = ["dp_ann_201904_an19_ymazefos_020719_pfc_z380-399_dp_ann.RoiSet.zip",
        "dp_ann_201904_an19_ymazefos_020719_cortex_z380-399_02_dp_ann.RoiSet.zip",
        "dp_ann_201904_an19_ymazefos_020719_thal_z350-369_dp_ann.RoiSet.zip",
        "dp_ann_an22_ymazecfos_z230-249_sm_cortex_dp_ann_RoiSet.zip",
        "dp_ann_201904_an19_ymazefos_020719_cortex_z350-369_dp_ann.RoiSet.zip",
        "dp_ann_201904_an12_ymazefos_020719_cortex_z371-390_dp_ann.RoiSet.zip",
        "tp_ann_201904_an4_ymazefos_020119_cortex_z200-219.RoiSet.zip",
        "dp_ann_201904_an22_ymazefos_020719_cb_z160-179_dp_ann.RoiSet.zip",
        "tp_ann_201904_an30_ymazefos_020719_striatum_z416-435.RoiSet.zip",
        "tp_ann_201904_an19_ymazefos_020719_pfc_z380-399.RoiSet.zip",
        "dp_ann_201904_an22_ymazefos_020719_midbrain_z150-169_dp_ann.RoiSet.zip",
        "tp_ann_201904_an10_ymzefos_020719_cortex_z280-279.RoiSet.zip",
        "dp_ann_201904_an12_ymazefos_020719_hypothal_z420-449_dp_ann.RoiSet.zip",
        "dp_ann_201904_an21_ymazefos_020719_hypothal_z450-469_dp_ann.RoiSet.zip",
        "dp_ann_an16_ymazecfos_z260-299_retrosplenial_dp_ann.Roiset.zip",
        "tp_ann_201904_an22_ymazefos_020719_pfc_z150-169.RoiSet.zip",
        "dp_ann_201904_an19_ymazefos_020719_cb_z380-399_dp_ann.RoiSet.zip",
        "tp_ann_201904_an4_ymazefos_020119_pfc_z200-219.RoiSet.zip"
        ]
roi_pths = [os.path.join(src, "raw_inputs/"+xx) for xx in rois]

dct = pickle.load(open(os.path.join(src, "clearmap_thresholds_sweep.p"), "rb"), encoding = "latin1")

#these will be zyx 
#note that importing it this way, the z dimension does not start from 0, but neither does clearmaps, so probably ok???
annotated_cells = np.array([np.array([[int(yy) for yy in xx[0].replace(".roi", "").split("-")] for xx in
                            read_roi_zip(roi_pth, include_roi_name=True)]) for roi_pth in roi_pths])

#voxel pair cutoff
cutoffs = [5,7]

for cutoff in cutoffs:
    #init dataframe
    print("cutoff: %s \n\n" % cutoff)
    df = pd.DataFrame()
    df["parameters"] = list(dct.keys())
    
    df["tp"] = np.zeros(len(df))
    df["fp"] = np.zeros(len(df))
    df["fn"] = np.zeros(len(df))
    
    df["f1"] = np.zeros(len(df))
    
    for n, fld in enumerate(list(dct.keys())):
        if n%100 == 0: print(n)
        
        detected_cells = np.asarray([v for k,v in dct[fld].items()])
        
        measures = [pairwise_distance_metrics(annotated_cells[i], detected_cells[i], cutoff = cutoff, verbose = False) for i
                    in range(len(detected_cells))] 
        
        tp = sum([measure[1] for measure in measures])
        fp = sum([measure[2] for measure in measures])
        fn = sum([measure[3] for measure in measures])
        
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
                
    df.to_csv(os.path.join(dst, "clearmap_performance_measures_voxelcutoff_%02d_v2.csv" % cutoff))
    
