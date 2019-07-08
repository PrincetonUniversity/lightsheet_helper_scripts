#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 16:10:57 2019

@author: wanglab
"""

import os, numpy as np
os.chdir("/jukebox/wang/zahra/python/lightsheet_py3")
from skimage.external import tifffile
import seaborn as sns, pandas as pd, matplotlib.pyplot as plt
import scipy, itertools
from skimage.exposure import equalize_hist, adjust_gamma
from tools.utils.io import listdirfull
from tools.analysis.network_analysis import make_structure_objects

sns.set_style('white')

#make inputs
src = "/jukebox/wang/Jess/lightsheet_output/201810_cfos/pooled_analysis/lobuleVI"
if not os.path.exists(src): os.mkdir(src)
flds = "/jukebox/wang/Jess/lightsheet_output/201810_cfos/processed"

#get files
lst = [fld for fld in listdirfull(flds) if os.path.exists(os.path.join(fld, 'Annotated_counts_60um_erosion.csv')) and fld[-7:-3] == "lob6" 
       or fld[-7:] == 'crus1_3' or fld[-7:] == 'crus1_1' or fld[-7:] == 'crus1_4' ]; lst.sort()

#conditions
nms = ['dadult_pc_crus1_1',
       'dadult_pc_crus1_3',
       'dadult_pc_crus1_4',
       'dadult_pc_lob6_13',
       'dadult_pc_lob6_14',
       'dadult_pc_lob6_15',
       'dadult_pc_lob6_16',
       'dadult_pc_lob6_17',
       'dadult_pc_lob6_18',
       'dadult_pc_lob6_19',
       'dadult_pc_lob6_20',
       'dadult_pc_lob6_21']

cond = ['Vector Control', 'Vector Control', 'Vector Control', 'DREADDs', 'Vector Control', 'DREADDs', 'Vector Control', 'DREADDs', 'DREADDs','DREADDs', 'DREADDs',
        'DREADDs']

conditions = {n:c for n,c in zip(nms, cond)}
pth = os.path.join(src, "cell_counts.csv")

df_pth = '/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx'
ann_pth = '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_60um_erosion.tif'
atl_pth = '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif'
#%%

#build structures class
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#%%
def generate_data_frame(conditions, lst, pth):
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
        df = pd.read_csv(fl+"/Annotated_counts_60um_erosion.csv")[1:] #remove previous headers
        print(nm, df.shape)
        df = df.replace(np.nan, "", regex=True)
        df["Brain"] = nm
        df["Condition"] = conditions[nm]
        bglst.append(df)
    
    df = pd.concat(bglst)
    df["counts"] = df["counts"].apply(int)

    #remove structures we don't care about!!!!!!!!!!!!!!!!!!!!!!
    sois = ["ventricular systems", "fiber tracts", "grooves"]
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        df = df[df.name != soi.name]
        progeny = [str(xx.name) for xx in soi.progeny]
        for progen in progeny:
            df = df[df.name != progen]
    
    df.drop(columns = ["Unnamed: 0"]).to_csv(pth, index = None)
    
    return pth

#run
csv_pth = generate_data_frame(conditions, lst, pth)

#%%
def generate_percent_counts_and_density_per_region(src, csv_pth):
    
    """ generates another column in the dataframe that is just # counts / total counts in brain """
    
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
    df["volume"] = df[df.voxels_in_structure > 0].apply(lambda x: x.voxels_in_structure*(scale**3), 1)
    df["density"] = df[df.voxels_in_structure > 0].apply(lambda x:x.counts/(float(x.voxels_in_structure*(scale**3))), 1)
    #export
    df.to_csv(os.path.join(src, "cell_counts_percents_density.csv"))
    
    return os.path.join(src, "cell_counts_percents_density.csv")

#run
percent_density_csv_pth = generate_percent_counts_and_density_per_region(src, csv_pth)
#%%
""" pools regions together based on allen name """    
 
#give list of structures you want to pool
pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/structures.csv"

sois = pd.read_csv(pth)
sois = [xx[0] for xx in sois.values]
   
#set variables
orgdf = pd.read_csv(percent_density_csv_pth).drop(columns = ["Unnamed: 0"])

#init
tdf = orgdf.copy()

#we have to do this per animal
#get the brains
brains = tdf.Brain.unique()

#make a big dicitionary to put it all together
pooled_counts_an = {}

#iterate through the conditions
for an in brains:
    df = tdf[tdf.Brain == an]
    condition = df["Condition"].unique()
    pooled_counts = {}
    
    #loop through and find downstream structures of a given parent
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print(soi.name)
        progeny = [str(xx.name) for xx in soi.progeny]
        counts = [] #store counts in this list
        val = df.loc[(df.name == soi.name), "percent"].values
        if val.shape[0] > 0: counts.append(val[0])
        if len(progeny) > 0:
            for progen in progeny:
                val = df.loc[(df.name == progen), "percent"].values
                if val.shape[0] > 0: counts.append(val[0])
                
        #sum counts in parent structures            
        pooled_counts[soi.name] = np.sum(np.asarray(counts))
    
    #fill other details and add to big dict
    pooled_counts_an[an] = pooled_counts
    pooled_counts["group"] = condition[0]
    
#make into giant dataframe
main_df = pd.DataFrame.from_dict(pooled_counts_an, orient = "index")
#sort by group so makes more sense
main_df = main_df.sort_values(by=["group"])
main_df.index.name = "animal"

main_df.to_csv(os.path.join(src, "select_structures_percent_counts_for_visualization.csv"))

rotate_df = pd.DataFrame()
struct = [list(itertools.repeat(xx, len(main_df))) for xx in main_df.columns.values[:-1]]
struct = pd.Series(list(itertools.chain.from_iterable(struct)))
rotate_df["name"] = struct

vals = [pd.Series(main_df[xx].values) for xx in main_df.columns.values[:-1]]    
rotate_df["percent"] = pd.concat(vals, ignore_index = True)

ans = list(itertools.repeat(main_df.index.values, len(struct)))
ans = pd.Series(list(itertools.chain.from_iterable(ans)))
rotate_df["animal"] = ans

groups = list(itertools.repeat(main_df.group.values, len(struct)))
groups = pd.Series(list(itertools.chain.from_iterable(groups)))
rotate_df["condition"] = groups

#save
rotate_df.to_csv(os.path.join(src, "select_structures_percent_counts_for_plots.csv"), index = False)

print("saved in :{}".format(src))

#%%
#anova for cell counts, percents, and density across all conditions, per structure
from scipy.stats import f_oneway

#do first for all structures
df = pd.read_csv(os.path.join(src, 'select_structures_percent_counts_for_plots.csv'))
#df = pd.read_csv(os.path.join(path,'select_structures_percent_counts_for_plots.csv'))

df_anova = pd.DataFrame()
df_anova["name"] = np.unique(df["name"].values)

for nm in np.unique(df.name.values)[:-1]: #only gets unique names

    f, pval = f_oneway(df[(df.name == nm) & (df.condition == "DREADDs")]["percent"].values, 
                 df[(df.name == nm) & (df.condition == "Vector Control")]["percent"].values)
    
    df_anova.loc[(df_anova["name"] == nm), "anova_percent_counts_f"] = f
    df_anova.loc[(df_anova["name"] == nm), "anova_percent_counts_pval"] = pval
        
df_anova.to_csv(os.path.join(src, "one_way_anova_select_structures.csv"))