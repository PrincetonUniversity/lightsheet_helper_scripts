#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:35:17 2020

@author: wanglab
"""


import os, pandas as pd, numpy as np, json, pickle

src = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"

data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/thal_model_data_contra_allen.p"
data = pickle.load(open(data_pth, "rb"), encoding = "latin1")

#set dest
dst = "/home/wanglab/Desktop"

#set the appropritate variables
primary_pool = data["primary_pool"]
ak_pool = data["ak_pool"]

cells_regions_pth = os.path.join(src, "thal_contra_counts_23_brains.csv")
cells_regions = pd.read_csv(cells_regions_pth)
brains = list(cells_regions.columns)[1:]

#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.020
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

def get_progeny(dic,parent_structure,progeny_list):
    """ 
    ---PURPOSE---
    Get a list of all progeny of a structure name.
    This is a recursive function which is why progeny_list is an
    argument and is not returned.
    ---INPUT---
    dic                  A dictionary representing the JSON file 
                         which contains the ontology of interest
    parent_structure     The structure
    progeny_list         The list to which this function will 
                         append the progeny structures. 
    """
    if 'msg' in list(dic.keys()): dic = dic['msg'][0]
    
    name = dic.get('name')
    children = dic.get('children')
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get('name')
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get('name')
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 

#get progeny of all large structures
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for all of neocortex
sois = ["Thalamus, sensory-motor cortex related", "Thalamus, polymodal association cortex related"]

#first calculate counts across entire nc region
contra_counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    contra_counts_per_struct.append(np.array(counts).sum(axis = 0))
contra_counts_per_struct = np.array(contra_counts_per_struct)


#%%

cells_regions_pth = os.path.join(src, "thal_ipsi_counts_23_brains.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])

#first calculate counts across entire nc region
ipsi_counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    ipsi_counts_per_struct.append(np.array(counts).sum(axis = 0))
ipsi_counts_per_struct = np.array(ipsi_counts_per_struct)

#%%
#get contra/ipsi ratios

ratio = contra_counts_per_struct.T/ipsi_counts_per_struct.T
mean_ratio = np.mean(ratio, axis = 0)
std_ratio = np.std(ratio, axis = 0)

#separate by injection, vermis vs. hemisphere

func = lambda xx: 0 if xx < 3 else 1
#prv
primary_pool_vh = np.array([func(xx) for xx in primary_pool])
ratio_vermis = ratio[np.where(primary_pool_vh == 0)]
ratio_hem = ratio[np.where(primary_pool_vh == 1)]

mean_ratio_vermis = np.mean(ratio_vermis, axis = 0)
mean_ratio_hem = np.mean(ratio_hem, axis = 0)
std_ratio_vermis = np.std(ratio_vermis, axis = 0)
std_ratio_hem = np.std(ratio_hem, axis = 0)