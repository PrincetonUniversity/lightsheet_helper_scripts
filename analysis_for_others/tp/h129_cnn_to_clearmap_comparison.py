#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:40:40 2019

@author: wanglab
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from tools.conv_net.utils.io import pairwise_distance_metrics, read_roi_zip
from tools.utils.io import listdirfull

pth = "/jukebox/wang/zahra/cnn_to_clearmap_comparison/roc_curve/"
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
    
#%%

dfs = listdirfull(dst)
dfs = [pd.read_csv(xx).sort_values(by = ["parameters"]) for xx in dfs if ".csv" in xx]
tps = np.array([df.tp.values for df in dfs]).sum(axis = 0)
fps = np.array([df.fp.values for df in dfs]).sum(axis = 0)
fns = np.array([df.fn.values for df in dfs]).sum(axis = 0)
precisions_cm = [(tp/(tp+fps[i])) for i, tp in enumerate(tps)]
#add 1 for plot
precisions_cm.append(1); precisions_cm = np.array(precisions_cm)
recalls_cm = [(tp/(tp+fns[i])) for i, tp in enumerate(tps)]
#add 0 for plot
recalls_cm.append(0); recalls_cm = np.array(recalls_cm)

src = "/jukebox/wang/zahra/conv_net/training/h129/experiment_dirs/20181115_zd_train/precision_recall_curve_295590.csv"
df = pd.read_csv(src)
precisions_cn = df["precision"].values
recalls_cn = df["recall"].values

#calculate
roc_auc_cm = np.trapz(precisions_cm, x = [1-xx for xx in recalls_cm])
roc_auc_cn = np.trapz(precisions_cn, x = [1-xx for xx in recalls_cn])

plt.figure()
plt.plot([1-xx for xx in recalls_cm], precisions_cm, color="darkorange", lw=1 , label="ClearMap ROC (area = %0.3f)" % roc_auc_cm)
plt.plot([1-xx for xx in recalls_cn], precisions_cn, color="darkred", lw=1 , label="CNN ROC (area = %0.3f)" % roc_auc_cn)
plt.plot([0, 1], [0, 1], color="navy", linestyle="--")
plt.xlim([0, 1])
plt.ylim([0, 1.05])
plt.xlabel("1 - Recall") 
plt.ylabel("Precision")
plt.title("Precision-Recall curve")
plt.legend(loc="lower right")
plt.savefig("/jukebox/wang/zahra/cnn_to_clearmap_comparison/roc_curve/results/cnn_to_clearmap_comparison_pr_curve.pdf")

