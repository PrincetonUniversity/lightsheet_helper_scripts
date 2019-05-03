#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:40:12 2019

@author: wanglab
"""

from __future__ import division
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

src = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection"
flds = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/processed"

#get files
lst = [fld for fld in listdirfull(flds) if os.path.exists(os.path.join(fld, "Annotated_counts.csv"))]

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
ann_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_75um_erosion_100um_ventricular_erosion.tif"
atl_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/average_template_25_sagittal_forDVscans_z_thru_240.tif"

#%%

#generate data frame
bglst=[]
for fl in lst:
    #extract out info
    nm = os.path.basename(fl)
    #make dataframe
    df = pd.read_csv(fl+"/Annotated_counts.csv" , error_bad_lines=False, names=["Index", "Count", "struct", "sub", "layer"])
    print(nm, df.shape)
    df = df.replace(np.nan, "", regex=True)
    f1 = df.struct.values
    f2 = df["sub"].values
    f3 = df.layer.values
    def remove(l, item):
        while item in l:
            l.remove(item)
        return l
    structures = [",".join(remove(list(l),"")) for l in zip(f1,f2,f3)]
    structures = [xx.replace("\xef\xbf\xbd", "e") for xx in structures]
    structures = [xx[1:] for xx in structures] #remove space from the beginning of each name this is important
    df["Structure"] = structures
    df["Brain"] = nm
    df["Condition"] = conditions[nm]
    #df["Stim_source"] = source[nm]
    del df["struct"]
    del df["sub"]
    del df["layer"]
    bglst.append(df)

df = pd.concat(bglst)
df["Count"] = df["Count"].apply(int)
df.to_csv("/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/cell_counts_dataframe_75um_erosion_100um_ventricular_erosion.csv", index = None)
#%%
#make percent counts per brain FOR BOTH COHORTS
csv_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/cell_counts_dataframe_75um_erosion_100um_ventricular_erosion.csv"
src = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection"
lut_pth = "/jukebox/wang/zahra/clearmap_cluster_copy/ClearMap/Data/ARA2_annotation_info_w_voxel_counts.csv"
ann_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans.nrrd"

#read csv
df = pd.read_csv(csv_pth, index_col = None).drop(columns = ["Index"])
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
df.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_w_percent_density.csv"), index = None)

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

#build structures class
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#set variables
csv_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_w_percent_density.csv"

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
#save both as csv and pickle for visualisation
tdf.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_with_progenitors_w_percents_density.csv"), header = True, index = None)

print("saved in :{}".format(src))


#%%
#############################################################CUSTOM SCRIPT#####################################################################################

#aggregate counts
df_new = tdf.groupby(["progenitor", "Brain"])["Count", "percent", "density"].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_projenitor_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"), header = True)

print("saved in :{}".format(src))

#%%

#set variables
df = pd.read_csv(csv_pth)
sois = ["Hippocampal formation", "Ammon"s horn", "Dentate gyrus", "Midbrain, sensory related", "Midbrain, motor related", 
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

#save as csv 
tdf.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_with_parents_w_percents_density.csv"), header = True, index = None)

print("saved in :{}".format(src))
    
#%%
#aggregate counts
df_new = tdf.groupby(["parent", "Brain"])["Count", "percent", "density"].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_parents_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"), header = True)

print("saved in :{}".format(src))

#%%

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

#save as csv 
tdf.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_with_children_w_percents_density.csv"), header = True, index = None)

print("saved in :{}".format(src))    
 
#%%
#aggregate counts
df_new = tdf.groupby(["child", "Brain"])["Count", "percent", "density"].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_children_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"), header = True)

print("saved in :{}".format(src))

#%%
#again...
pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/cell_counts_dataframe_with_children_w_percents_density.csv"
tdf = pd.read_csv(pth, index_col = None).drop(columns = ["Unnamed: 0"])
#aggregate counts
df_new = tdf.groupby(["child", "Brain"])["Count", "percent", "density"].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_children_cell_counts_dataframe_w_percents_density.csv"), header = True)

print("saved in :{}".format(src))
#%%
#concantenate all dfs
df1 = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/summed_children_cell_counts_dataframe_w_percents_density.csv"
df2 = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/summed_parents_cell_counts_dataframe_w_percents_density.csv"
df3 = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/summed_projenitor_cell_counts_dataframe_w_percents_density.csv"

#formatting...
df1 = pd.read_csv(df1, index_col = None)
df1["structure"] = df1.child
df1 = df1.drop(columns = ["child"])
df2 = pd.read_csv(df2, index_col = None)
df2["structure"] = df2.parent
df2 = df2.drop(columns = ["parent"])
df3 = pd.read_csv(df3, index_col = None)
df3["structure"] = df3.progenitor
df3 = df3.drop(columns = ["progenitor"])

df_no_erosion = pd.concat([df1, df2, df3])

#%%
df4 = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/summed_children_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"
df5 = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/summed_parents_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"
df6 = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/summed_projenitor_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"

#formatting...
df4 = pd.read_csv(df4, index_col = None)
df4["structure"] = df4.child
df4 = df4.drop(columns = ["child"])

df5 = pd.read_csv(df5, index_col = None)
df5["structure"] = df5.parent
df5 = df5.drop(columns = ["parent"])

df6 = pd.read_csv(df6, index_col = None)
df6["structure"] = df6.progenitor
df6 = df6.drop(columns = ["progenitor"])

df_erosion = pd.concat([df4, df5, df6])

#%%

df = df_no_erosion.copy()

df["erosion_Count"] = df_erosion.Count
df["erosion_percent"] = df_erosion.percent
df["erosion_density"] = df_erosion.density

df = df[["structure", "Brain", "Count", "percent", "density",
       "erosion_Count", "erosion_percent", "erosion_density"]]

df.to_csv("/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/analysis/cell_detection/composite_cell_detection_results_6mo_cohort_w_percent_density_erosion.csv", index = None)


#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#%%

src = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection"
flds = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/processed"

#145 none is left side of brain?

#get files
lst = [fld for fld in listdirfull(flds) if os.path.exists(os.path.join(fld, "Annotated_counts.csv"))]

#conditions
nms = ["none_c6_488_64","pbs_c6_488_647","v143_none_c1_1","v143_rt_c1_2_4","v144_none_c4_1","v144_rt_c4_2_4",
       'v145_none_c3_1', "v145_rt_c3_2_4","v190_none_c5_1","v190_rt_c5_2_4","v75_none_c2_1_","v75_rt_c2_2_48"]
cond = ["Control", "Control", "v143", "v143", "v144", "v144", "v145", "v145", "v190", "v190", "v75", "v75"]
conditions = {n:c for n,c in zip(nms, cond)}

df_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/ARA2_annotation_info_w_voxel_counts.csv"
atl_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/average_template_25_sagittal_forDVscans_z_thru_240.tif"

#%%

#generate data frame
bglst=[]
for fl in lst:
    #extract out info
    nm = os.path.basename(fl)
    #make dataframe
    df = pd.read_csv(fl+"/Annotated_counts.csv" , error_bad_lines=False, names=["Index", "Count", "struct", "sub", "layer"])
    print(nm, df.shape)
    df = df.replace(np.nan, "", regex=True)
    f1 = df.struct.values
    f2 = df["sub"].values
    f3 = df.layer.values
    def remove(l, item):
        while item in l:
            l.remove(item)
        return l
    structures = [",".join(remove(list(l),"")) for l in zip(f1,f2,f3)]
    structures = [xx.replace("\xef\xbf\xbd", "e") for xx in structures]
    structures = [xx[1:] for xx in structures] #remove space from the beginning of each name this is important
    df["Structure"] = structures
    df["Brain"] = nm
    df["Condition"] = conditions[nm]
    #df["Stim_source"] = source[nm]
    del df["struct"]
    del df["sub"]
    del df["layer"]
    bglst.append(df)

df = pd.concat(bglst)
df["Count"] = df["Count"].apply(int)
df.to_csv("/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/cell_counts_dataframe_75um_erosion_100um_ventricular_erosion.csv", index = None)

#%%
#make percent counts per brain FOR BOTH COHORTS
csv_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/cell_counts_dataframe_75um_erosion_100um_ventricular_erosion.csv"
src = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection"
lut_pth = "/jukebox/wang/zahra/clearmap_cluster_copy/ClearMap/Data/ARA2_annotation_info_w_voxel_counts.csv"
ann_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_75um_erosion_100um_ventricular_erosion.tif"

#read csv
df = pd.read_csv(csv_pth, index_col = None).drop(columns = ["Index"])
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
df.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_w_percent_density.csv"), index = None)

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

#build structures class
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#set variables
csv_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_w_percent_density.csv"

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
#save both as csv and pickle for visualisation
tdf.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_with_progenitors_w_percents_density.csv"), header = True, index = None)

print("saved in :{}".format(src))


#%%
#aggregate counts
df_new = tdf.groupby(["progenitor", "Brain"])["Count", "percent", "density"].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_progenitor_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"), header = True)

print("saved in :{}".format(src))

#%%

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

#save as csv 
tdf.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_with_parents_w_percents_density.csv"), header = True, index = None)

print("saved in :{}".format(src))
    
#%%
#aggregate counts
df_new = tdf.groupby(["parent", "Brain"])["Count", "percent", "density"].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_parents_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"), header = True)

print("saved in :{}".format(src))

#%%

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

#save as csv 
tdf.to_csv(os.path.join(src, "cell_counts_dataframe_75um_erosion_100um_ventricular_erosion_with_children_w_percents_density.csv"), header = True, index = None)

print("saved in :{}".format(src))    
 
#%%
#aggregate counts
df_new = tdf.groupby(["child", "Brain"])["Count", "percent", "density"].sum()

#save as csv 
df_new.to_csv(os.path.join(src, "summed_children_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"), header = True)

print("saved in :{}".format(src))

#%%
#concantenate all dfs
df1 = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/summed_children_cell_counts_dataframe_w_percents_density.csv"
df2 = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/summed_parents_cell_counts_dataframe_w_percents_density.csv"
df3 = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/summed_progenitor_cell_counts_dataframe_w_percents_density.csv"

#formatting...
df1 = pd.read_csv(df1, index_col = None)
df1["structure"] = df1.child
df1 = df1.drop(columns = ["child"])
df2 = pd.read_csv(df2, index_col = None)
df2["structure"] = df2.parent
df2 = df2.drop(columns = ["parent"])
df3 = pd.read_csv(df3, index_col = None)
df3["structure"] = df3.progenitor
df3 = df3.drop(columns = ["progenitor"])

df_no_erosion = pd.concat([df1, df2, df3])

#%%
df4 = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/summed_children_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"
df5 = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/summed_parents_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"
df6 = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/summed_progenitor_cell_counts_dataframe_w_percents_density_75um_erosion_100um_ventricular_erosion.csv"

#formatting...
df4 = pd.read_csv(df4, index_col = None)
df4["structure"] = df4.child
df4 = df4.drop(columns = ["child"])

df5 = pd.read_csv(df5, index_col = None)
df5["structure"] = df5.parent
df5 = df5.drop(columns = ["parent"])

df6 = pd.read_csv(df6, index_col = None)
df6["structure"] = df6.progenitor
df6 = df6.drop(columns = ["progenitor"])

df_erosion = pd.concat([df4, df5, df6])

#%%

df = df_no_erosion.copy()

df["erosion_Count"] = df_erosion.Count
df["erosion_percent"] = df_erosion.percent
df["erosion_density"] = df_erosion.density

df = df[["structure", "Brain", "Count", "percent", "density",
       "erosion_Count", "erosion_percent", "erosion_density"]]

df.to_csv("/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/analysis/cell_detection/composite_cell_detection_results_1mo_cohort_w_percent_density_erosion.csv", index = None)
