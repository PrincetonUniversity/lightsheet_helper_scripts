#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 17:19:59 2020

@author: wanglab
"""

import tifffile, os, numpy as np, pandas as pd

if __name__ == "__main__":
    
    #read transformed file
    src = "/jukebox/LightSheetData/kocher-bee/volume_analysis/template_to_brain"
    trsfmpths = [os.path.join(src,xx,"result.tif") for xx in os.listdir(src)]
    annpth = "/jukebox/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment_croppedZ.tif"
    ann = tifffile.imread(annpth)
    #make dataframe to store values
    df = pd.DataFrame()
    #find all unique annotated regions
    ids = np.unique(ann)
    #iterate through all ids and find total no. of voxels in each region of the atlas
    df["template"] = np.zeros(len(ids))
    for iid in ids:
        print(iid)
        voxels = np.sum((ann==iid).astype(int))
        df.loc[iid, "template"] = voxels
        
    #iterate through all registered volumes
    for trsfmpth in trsfmpths:
        trsfm = tifffile.imread(trsfmpth)
        brain = os.path.basename(os.path.dirname(trsfmpth))
        df[brain] = np.zeros(len(ids))
        print(brain)
        for iid in ids:
            #iterate through all annotated regions and find # of voxels per region
            print(iid)
            voxels = np.sum((trsfm==float(iid)).astype(int))
            df.loc[iid, brain] = voxels
        
    #export
    df.to_csv(os.path.join(src, "voxels_per_region.csv"))
