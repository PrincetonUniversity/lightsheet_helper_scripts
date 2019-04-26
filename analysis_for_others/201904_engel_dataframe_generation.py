#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:40:12 2019

@author: wanglab
"""

import os, numpy as np
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from skimage.external import tifffile
import seaborn as sns, pandas as pd, matplotlib.pyplot as plt
import scipy, itertools
from skimage.exposure import equalize_hist, adjust_gamma
from tools.utils.io import listdirfull
from tools.analysis.network_analysis import make_structure_objects
from tools.utils.overlay import tile
os.chdir("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.cluster.utils import load_kwargs
import ClearMap.IO.IO as io
import ClearMap.Analysis.Statistics as stat
import ClearMap.Alignment.Resampling as rsp

sns.set_style("white")

#%%
src = "/home/wanglab/mounts/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection"
flds = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/processed"

#get files
lst = [fld for fld in listdirfull(flds) if os.path.exists(os.path.join(fld, "Annotated_counts_75um_erosion.csv"))]

#conditions
nms = ["buffer",
     "control",
     "v143_3",
     "v143_4",
     "v144_3",
     "v144_4",
     "v145_3",
     "v145_4",
     "v190_3",
     "v190_4",
     "v75_3",
     "v75_4"]

pths = [os.path.join(flds, nm) for nm in nms]

cond = ["Buffer", "Control", "v143", "v143", "v144", "v144", "v145", "v145", "v190", "v190", "v75", "v75"]
conditions = {n:c for n,c in zip(nms, cond)}
path_per_condition = {n:c for n,c in zip(pths, cond)}

df_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/ARA2_annotation_info_w_voxel_counts.csv"
ann_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans_z_thru_240.nrrd"
atl_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/average_template_25_sagittal_forDVscans_z_thru_240.tif"

#%%

#redo cell counts dataframe
#generate data frame
bglst=[]
for fl in lst:
    #extract out info
    nm = os.path.basename(fl)
    #make dataframe
    df = pd.read_csv(fl+'/Annotated_counts_75um_erosion.csv' , error_bad_lines = False)
    print(nm, df.shape)
    df = df.replace(np.nan, '', regex=True)
    f1 = df.struct.values
    f2 = df['sub'].values
    f3 = df.layer.values
    def remove(l, item):
        while item in l:
            l.remove(item)
        return l
    structures = [','.join(remove(list(l),'')) for l in zip(f1,f2,f3)]
    structures = [xx.replace('\xef\xbf\xbd', 'e') for xx in structures]
    structures = [xx[1:] for xx in structures] #remove space from the beginning of each name this is important
    df['Structure'] = structures
    df['Brain'] = nm
    df['Condition'] = conditions[nm]
    #df['Stim_source'] = source[nm]
    del df['struct']
    del df['sub']
    del df['layer']
    bglst.append(df)

df = pd.concat(bglst)
df['Count'] = df['Count'].apply(int)
df.to_csv('/home/wanglab/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/cell_counts_dataframe.csv')


#%%
#make percent counts per brain FOR BOTH COHORTS
csv_pth = '/home/wanglab/mounts/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/cell_counts_dataframe.csv'
src = "/home/wanglab/mounts/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection"
lut_pth = "/home/wanglab/mounts/wang/zahra/clearmap_cluster_copy/ClearMap/Data/ARA2_annotation_info_w_voxel_counts.csv"
ann_pth = '/home/wanglab/mounts/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans.nrrd'

#read csv
df = pd.read_csv(csv_pth)
lut = pd.read_csv(lut_pth)

#make voxel counts column
structs = df.Structure.values

#initialise
df["voxels_in_structure"] = 0

for nm in structs:
    if not len(lut.loc[lut.name == nm, "voxels_in_structure"].values) == 0:
        df.loc[df.Structure == nm, "voxels_in_structure"] = lut.loc[lut.name == nm, "voxels_in_structure"].values[0]

#set resolution of atlas
scale = 0.025 ##mm/voxel

#get all brain names
brains = np.unique(df.Brain.values)

#save total counts in dict
total_counts = {}
percents = {}

#for each brain, get total counts
for brain in brains:
    total_counts[brain] = df[df.Brain == brain].Count.sum(0)    
       
percents = [df[df.Brain == brain].Count.apply(lambda x: (x/total_counts[brain])*100).astype("float64") for brain in brains]

#concantenate together
df_percent = pd.concat(percents)

df["percent"] = df_percent
df["volume"] = df[df.voxels_in_structure > 0].apply(lambda x: x.voxels_in_structure*(scale**3), 1)
df["density"] = df[df.voxels_in_structure > 0].apply(lambda x:x.Count/(float(x.voxels_in_structure*(scale**3))), 1)
#export
df.to_csv(os.path.join(src, "cell_counts_dataframe_w_percents_density.csv"))

#%%
#helper functions
def correct_cm_to_sturctures(struc):
    """function to correct naming issues
    """
    if struc == "Anterior cingulate area, ventral part, layer 6a": struc = "Anterior cingulate area, ventral part, 6a"
    if struc == "Anterior cingulate area, ventral part, layer 6b": struc = "Anterior cingulate area, ventral part, 6b"
    if struc == "Simple lobule": struc = "Simplex lobule"
    if struc == "Primary somatosensory area, barrel field, layer 5 ": struc = "Primary somatosensory area, barrel field, layer 5"
    return struc

def correct_sturctures_to_cm(struc):
    """function to correct naming issues
    """
    if struc == "Anterior cingulate area, ventral part, 6a": struc = "Anterior cingulate area, ventral part, layer 6a"
    if struc == "Anterior cingulate area, ventral part, 6b": struc = "Anterior cingulate area, ventral part, layer 6b"
    if struc == "Simplex lobule": struc = "Simple lobule"
    if struc == "Primary somatosensory area, barrel field, layer 5": struc = "Primary somatosensory area, barrel field, layer 5 "
    return struc  

#################################################################################################################################################################


def generate_progenitors_normalised_structures_list(df_pth, ann_pth, csv_pth):
    """
    generates counts normalised by volume, correct some errors in look up table
    #TODO (zahra): ask tom if this is necessary for PMA
    """
    #build structures class
    structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
    
    #set variables
    df = pd.read_csv(csv_pth)
    sois = ["Cerebral cortex", "Thalamus", "Hypothalamus", "Midbrain", "Hindbrain", 
            "Cerebellum", "Striatum", "Pallidum", "fiber tracts", "Olfactory areas"
            ]

    #add in progenitor column and add in volumes for area cell counts
    vols = pd.read_csv("/jukebox/wang/zahra/clearmap_cluster_copy/ClearMap/Data/ARA2_annotation_info_w_voxel_counts.csv")[["voxels_in_structure", "name"]]
    tdf = df.copy()
    tdf["progenitor"] = "empty"
    tdf["Volume"] = 0.0
    scale_factor = .025 #mm/voxel
    
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print(soi.name)
        progeny = [str(xx.name) for xx in soi.progeny]
        if len(progeny) > 0:
            for progen in progeny:
                print(progen)
                progen = correct_sturctures_to_cm(progen)
                tdf.loc[tdf["Structure"] == progen,"progenitor"]=soi.name
                if len(vols[vols["name"] == progen]["voxels_in_structure"])>0:
                    tdf.loc[tdf["Structure"] == progen,"Volume"] = vols[vols["name"] == progen]["voxels_in_structure"].values[0]*scale_factor
        else:
            tdf.loc[tdf["Structure"] == soi.name,"progenitor"] = soi.name
            if len(vols[vols["name"] == soi.name]["voxels_in_structure"])>0:
                tdf.loc[tdf["Structure"] == soi.name,"Volume"] = vols[vols["name"] == soi.name]["voxels_in_structure"].values[0]*scale_factor
            
    
    #drop non progen
    tdf = tdf[(tdf["progenitor"] != "empty") & (tdf["Volume"] != 0.0)]
    
    #add normalised column
    tdf["count_normalized_by_volume"] = tdf.apply(lambda x:x["Count"]/float(x["Volume"]), 1)
    
    #drop unnamed
    tdf = tdf.drop(columns = ["Index"])
    #save both as csv and pickle for visualisation
    tdf.to_csv(os.path.join(src, "cell_counts_dataframe_with_progenitors_w_percents_density.csv"), header = True, index = None)
    
    print("saved in :{}".format(src))
    
    return tdf

#run
tdf = generate_progenitors_normalised_structures_list(df_pth, ann_pth, csv_pth)

#%%
#############################################################CUSTOM SCRIPT#####################################################################################

#aggregate counts
df_new = tdf.groupby(["progenitor", "Brain"])["Count", 'percent', 'density'].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_projenitor_cell_counts_dataframe_w_percents_density.csv"), header = True)

print("saved in :{}".format(src))

#%%
def generate_parents_normalised_structures_list(df_pth, ann_pth, csv_pth):
    """
    generates counts normalised by volume, correct some errors in look up table
    #TODO (zahra): ask tom if this is necessary for PMA
    """
    #build structures class
    structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
    
    #set variables
    df = pd.read_csv(csv_pth)
    sois = ["Hippocampal formation", "Ammon's horn", "Dentate gyrus", "Midbrain, sensory related", "Midbrain, motor related", 
            "Midbrain, behavioral state related", "Cerebellar cortex", "Cerebellar nuclei", "Primary motor area",
            "Secondary motor area", "Primary somatosensory area", "Supplemental somatosensory area", "Visual areas"]

    #add in progenitor column and add in volumes for area cell counts
    vols = pd.read_csv("/jukebox/wang/zahra/clearmap_cluster_copy/ClearMap/Data/ARA2_annotation_info_w_voxel_counts.csv")[["voxels_in_structure", "name"]]
    tdf = df.copy()
    tdf["parent"] = "empty"
    tdf["Volume"] = 0.0
    scale_factor = .025 #mm/voxel
    
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print(soi.name)
        progeny = [str(xx.name) for xx in soi.progeny]
        if len(progeny) > 0:
            for progen in progeny:
                progen = correct_sturctures_to_cm(progen)
                tdf.loc[tdf["Structure"] == progen,"parent"]=soi.name
                if len(vols[vols["name"] == progen]["voxels_in_structure"])>0:
                    tdf.loc[tdf["Structure"] == progen,"Volume"] = vols[vols["name"] == progen]["voxels_in_structure"].values[0]*scale_factor
            
    
    #drop non progen
    tdf = tdf[(tdf["parent"] != "empty") & (tdf["Volume"] != 0.0)]
    
    #add normalised column
    tdf["count_normalized_by_volume"] = tdf.apply(lambda x:x["Count"]/float(x["Volume"]), 1)
    
    #drop unnamed
    tdf = tdf.drop(columns = ["Index"])
    #save as csv 
    tdf.to_csv(os.path.join(src, "cell_counts_dataframe_with_parents_w_percents_density.csv"), header = True, index = None)
    
    print("saved in :{}".format(src))
    
    return tdf

#run
tdf = generate_parents_normalised_structures_list(df_pth, ann_pth, csv_pth)

#%%
#############################################################CUSTOM SCRIPT#####################################################################################

#aggregate counts
df_new = tdf.groupby(["parent", "Brain"])["Count", 'percent', 'density'].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_parents_cell_counts_dataframe_w_percents_density.csv"), header = True)

print("saved in :{}".format(src))

#%%
def generate_children_normalised_structures_list(df_pth, ann_pth, csv_pth):
    """
    generates counts normalised by volume, correct some errors in look up table
    #TODO (zahra): ask tom if this is necessary for PMA
    """
    #build structures class
    structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
    
    #set variables
    df = pd.read_csv(csv_pth)
    sois = ["Caudoputamen", "corpus callosum", "Nucleus accumbens"]

    #add in progenitor column and add in volumes for area cell counts
    vols = pd.read_csv("/jukebox/wang/zahra/clearmap_cluster_copy/ClearMap/Data/ARA2_annotation_info_w_voxel_counts.csv")[["voxels_in_structure", "name"]]
    tdf = df.copy()
    tdf["child"] = "empty"
    tdf["Volume"] = 0.0
    scale_factor = .025 #mm/voxel
    
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print(soi.name)
        tdf.loc[tdf["Structure"] == soi.name, "child"] = soi.name
        if len(vols[vols["name"] == soi.name]["voxels_in_structure"])>0:
            tdf.loc[tdf["Structure"] == soi.name,"Volume"] = vols[vols["name"] == soi.name]["voxels_in_structure"].values[0]*scale_factor
            
    
    #drop non progen
    tdf = tdf[(tdf["child"] != "empty")]
    
    #drop unnamed
    tdf = tdf.drop(columns = ["Index"])
    #save as csv 
    tdf.to_csv(os.path.join(src, "cell_counts_dataframe_with_children_w_percents_density.csv"), header = True, index = None)
    
    print("saved in :{}".format(src))
    
    return tdf

#run
tdf = generate_children_normalised_structures_list(df_pth, ann_pth, csv_pth)
#%%
#############################################################CUSTOM SCRIPT#####################################################################################

#aggregate counts
df_new = tdf.groupby(["child", "Brain"])['Count'].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_children_cell_counts_dataframe.csv"), header = True)

print("saved in :{}".format(src))
