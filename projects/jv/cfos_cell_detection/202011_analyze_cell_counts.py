#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:19:42 2020

@author: wanglab
"""

import os, numpy as np, pandas as pd, scipy, itertools, sys, json

def generate_data_frame(conditions, lst, pth, flnm):
    """ 
    used to make a pooled csv file of all cell counts in an experiment
    inputs:
        conditions: zip of file names + condition
        lst: list of file names run through analysis
        pth: path to save csv output
    """
    #generate data frame
    bglst=[]
    for fl in lst:
        #extract out info
        nm = os.path.basename(fl)
        #make dataframe
        df = pd.read_csv(fl+"/%s.csv" % flnm)[1:] #remove previous headers
        print(nm, df.shape)
        df = df.replace(np.nan, "", regex=True)
        df["Brain"] = nm
        df["Condition"] = conditions[nm]
        bglst.append(df)
    df = pd.concat(bglst)
    df["counts"] = df["counts"].apply(int)
    df.drop(columns = ["Unnamed: 0"]).to_csv(pth, index = None)
    
    return pth

def generate_percent_counts_and_density_per_region(src, csv_pth):
    
    """ generates another column in the dataframe that is just # counts / total counts in brain and density per region"""
    
    #read csv
    df = pd.read_csv(csv_pth)
    #set resolution of atlas
    scale = 0.020 ##mm/voxel
    #get all brain names
    brains = np.unique(df.Brain.values)
    #save total counts in dict
    total_counts = {}
    percents = {}
    #for each brain, get total counts
    for brain in brains:
        total_counts[brain] = df[df.Brain == brain].counts.sum(0)    
           
    percents = [df[df.Brain == brain].counts.apply(lambda x: (x/total_counts[brain])*100).astype("float64") for brain in brains]
    #concantenate together
    df_percent = pd.concat(percents)
    df["percent"] = df_percent
    df["density"] = df[df.voxels_in_structure > 0].apply(lambda x:x.counts/(float(x.voxels_in_structure*(scale**3))), 1)
    #export
    df.to_csv(os.path.join(src, "cell_counts.csv"))
    
    return os.path.join(src, "cell_counts.csv")

#%%
if __name__ == "__main__":
    
    #init source paths
    src = "/jukebox/wang/Jess/lightsheet_output/202010_cfos/pooled_analysis"
    if not os.path.exists(src): os.mkdir(src)
    flds = "/jukebox/wang/Jess/lightsheet_output/202010_cfos/processed"
    
    #get files
    lst = [os.path.join(flds, fld) for fld in os.listdir(flds) if os.path.exists(os.path.join(os.path.join(flds, fld), 
                                                         "Annotated_counts_60um_edge_80um_vent_erosion.csv"))]; lst.sort()
    #names and conditions
    nms = os.listdir(flds); nms.sort()
    cond = ["habituation","habituation","habituation","habituation","habituation","habituation","habituation","habituation",
            "habituation","habituation","acquisition_day1","acquisition_day1","acquisition_day1","acquisition_day1",
            "acquisition_day1","acquisition_day1","acquisition_day1","acquisition_day1","acquisition_day1","acquisition_day1",
            "acquisition_day2","acquisition_day2","acquisition_day2","acquisition_day2","acquisition_day2","acquisition_day2",
            "acquisition_day2","acquisition_day2","acquisition_day2","acquisition_day2"]
    conditions = {n:c for n,c in zip(nms, cond)}
    
    #set appropriate paths
    pth = os.path.join(src, "cell_counts.csv")
    df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
    ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
    atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
    
    #run
    csv_pth = generate_data_frame(conditions, lst, pth, "Annotated_counts_60um_edge_80um_vent_erosion")
    
    #run
    percent_density_csv_pth = generate_percent_counts_and_density_per_region(src, csv_pth)