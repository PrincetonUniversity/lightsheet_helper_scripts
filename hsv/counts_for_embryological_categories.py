#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:20:07 2020

@author: wanglab
"""

import os, pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, json, matplotlib as mpl, pickle as pckl

#get counts from gross embryological categories

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

def get_progeny(dic,parent_structure,progeny_list):
   
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

#%%

#start with thalamus

src = "/jukebox/wang/zahra/h129_contra_vs_ipsi"

data_pth = os.path.join(src, "data/thal_hsv_maps_contra_allen.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

primary_pool = data["primary_pool"]
frac_of_inj_pool = data["frac_of_inj_pool"]
ak_pool = data["ak_pool"]

df_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
cells_regions_pth_contra = os.path.join(src, "data/thal_contra_counts_23_brains_80um_ventric_erosion.csv")
cells_regions_pth_ipsi = os.path.join(src, "data/thal_ipsi_counts_23_brains_80um_ventric_erosion.csv")
dst = "/home/wanglab/Desktop"

cells_regions_contra = pd.read_csv(cells_regions_pth_contra)
cells_regions_ipsi = pd.read_csv(cells_regions_pth_ipsi)
cells_regions_ipsi = cells_regions_ipsi[cells_regions_contra.columns]
cells_regions = cells_regions_ipsi + cells_regions_contra #bilateral

#rename structure column
cells_regions["Structure"] = cells_regions_contra["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])

scale_factor = 0.025
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

#from disynaptic curation
brains = [#"20170410_tp_bl6_lob6a_ml_repro_01",
         "20160823_tp_bl6_cri_500r_02",
         "20180417_jg59_bl6_cri_03",
        # "20170207_db_bl6_crii_1300r_02",
        # "20160622_db_bl6_unk_01",
        # "20161205_tp_bl6_sim_750r_03",
         "20180410_jg51_bl6_lob6b_04",
        # "20170419_db_bl6_cri_rpv_53hr",
         "20170116_tp_bl6_lob6b_lpv_07",
        # "20170411_db_bl6_crii_mid_53hr",
         "20160822_tp_bl6_crii_1500r_06",
         "20160920_tp_bl6_lob7_500r_03",
        # "20170207_db_bl6_crii_rpv_01",
        # "20161205_tp_bl6_sim_250r_02",
        # "20161207_db_bl6_lob6a_500r_53hr",
        # "20170130_tp_bl6_sim_rlat_05",
         "20170115_tp_bl6_lob6b_500r_05",
       #  "20170419_db_bl6_cri_mid_53hr",
         "20161207_db_bl6_lob6a_850r_53hr",
       #  "20160622_db_bl6_crii_52hr_01",
         "20161207_db_bl6_lob6a_50rml_53d5hr",
       #  "20161205_tp_bl6_lob45_1000r_01",
         "20160801_db_l7_cri_01_mid_64hr"]


#get progeny of all large structures
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

sois = ["Cerebrum", #telencephalon,
        "Interbrain", #diencephalon
        "Midbrain", #mesencephalon
        "Pons", #metencephalon
        "Medulla"] #myelencephalon

#first calculate counts across entire region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    except:
        print(soi)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
    
counts_per_struct = np.array(counts_per_struct)

#voxels
vol = []
for soi in sois:
    progeny = []; counts = []; iids = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0])
    except:
        print(soi)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T)

#boxplots of counts
# plt.figure(figsize = (5,4))
# df = pd.DataFrame(counts_per_struct.T)
# df.columns = sois
# g = sns.stripplot(data = df,  color = "#99C7E0", orient = "h")
# sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
#             boxprops={'facecolor':'None'})

# plt.xlabel("Neurons")
# plt.ylabel("Region")
# plt.title("Thalamic timepoint")
# plt.savefig(os.path.join(dst, "gross_region_counts_boxplots_thal_tp.pdf"), bbox_inches = "tight")
# plt.close()

#boxplots of density
plt.figure(figsize = (5,4))
df = pd.DataFrame(density)
df.columns = sois
sns.stripplot(data = df,  color = "#99C7E0", orient = "h")
g = sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'})
g.set_xscale("log")

plt.xlabel("Neurons / mm$^3$")
plt.ylabel("Region")
plt.title("Thalamic timepoint")
plt.savefig(os.path.join(dst, "gross_region_density_boxplots_thal_tp.pdf"), bbox_inches = "tight")
plt.close()

#%%

#detailed structures
sois = ["Cerebral cortex", #telencephalon
        "Striatum",
        "Pallidum",
        "Thalamus",
        "Hypothalamus",
        "Ventral tegmental area", #diencephalon
        "Substantia nigra, compact part",
        "Substantia nigra, reticular part",
        "Red nucleus", #mesencephalon
        "Pontine gray", #metencephalon
        "Pontine reticular nucleus",
        "Inferior olivary complex"] #myelencephalon

#first calculate counts across entire region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    except:
        print(soi)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
    
counts_per_struct = np.array(counts_per_struct)

#voxels
vol = []
for soi in sois:
    progeny = []; counts = []; iids = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0])
    except:
        print(soi)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T)

# #boxplots of counts
# plt.figure(figsize = (5,5))
# df = pd.DataFrame(counts_per_struct.T)
# df.columns = sois
# g = sns.stripplot(data = df,  color = "dimgrey", orient = "h")
# sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
#             boxprops={'facecolor':'None'})

# # g.set(xscale = "log")
# plt.xlabel("Neurons")
# plt.ylabel("Region")
# plt.title("Thalamic timepoint")
# plt.savefig(os.path.join(dst, "detailed_region_counts_boxplots_thal_tp.pdf"), bbox_inches = "tight")
# plt.close()

# #boxplots of counts
# plt.figure(figsize = (5,5))
# g = sns.stripplot(data = df,  color = "dimgrey", orient = "h")
# sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
#             boxprops={'facecolor':'None'})

# g.set_xlim([-100, 2500])
# # g.set(xscale = "log")
# plt.xlabel("Neurons")
# plt.ylabel("Region")
# plt.title("Thalamic timepoint")
# plt.savefig(os.path.join(dst, "detailed_region_counts_boxplots_thal_tp_zoom.pdf"), bbox_inches = "tight")
# plt.close()

#boxplots of density
plt.figure(figsize = (5,5))
df = pd.DataFrame(density)
df.columns = sois
sns.stripplot(data = df,  color = "#99C7E0", orient = "h")
g = sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'})

g.set_xlim([-10, 220])
plt.xlabel("Neurons / mm$^3$")
plt.ylabel("Region")
plt.title("Thalamic timepoint")
plt.savefig(os.path.join(dst, "detailed_region_density_boxplots_thal_tp_zoom.pdf"), bbox_inches = "tight")
plt.close()
#%%

#neocortex

src = "/jukebox/wang/zahra/h129_contra_vs_ipsi"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx" #this is pma
cells_regions_pth_contra = os.path.join(src, "data/nc_contra_counts_33_brains_pma.csv")
cells_regions_pth_ipsi = os.path.join(src, "data/nc_ipsi_counts_33_brains_pma.csv")
dst = "/home/wanglab/Desktop"

cells_regions_contra = pd.read_csv(cells_regions_pth_contra)
cells_regions_ipsi = pd.read_csv(cells_regions_pth_ipsi)
cells_regions_ipsi = cells_regions_ipsi[cells_regions_contra.columns]
cells_regions = cells_regions_ipsi + cells_regions_contra #bilateral

#rename structure column
cells_regions["Structure"] = cells_regions_contra["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])

scale_factor = 0.025
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

brains = cells_regions.columns[:-1]

#get progeny of all large structures
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

sois = ["Cerebrum", #telencephalon,
        "Interbrain", #diencephalon
        "Midbrain", #mesencephalon
        "Pons", #metencephalon
        "Medulla"] #myelencephalon

#first calculate counts across entire region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    except:
        print(soi)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
    
counts_per_struct = np.array(counts_per_struct)

#voxels
vol = []
for soi in sois:
    progeny = []; counts = []; iids = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0])
    except:
        print(soi)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T) #remove thalamus

# #boxplots of counts
# plt.figure(figsize = (5,4))
# df = pd.DataFrame(counts_per_struct.T)
# df.columns = sois
# g = sns.stripplot(data = df,  color = "dimgrey", orient = "h")
# sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
#             boxprops={'facecolor':'None'})

# plt.xlabel("Neurons")
# plt.ylabel("Region")
# plt.title("Neocortical timepoint")
# plt.savefig(os.path.join(dst, "gross_region_counts_boxplots_nc_tp.pdf"), bbox_inches = "tight")
# plt.close()

#boxplots of density
plt.figure(figsize = (5,4))
df = pd.DataFrame(density)
df.columns = sois
g = sns.stripplot(data = df,  color = "#99C7E0", orient = "h")
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'})

plt.xlabel("Neurons/ mm$^3$")
plt.ylabel("Region")
plt.title("Neocortical timepoint")
plt.savefig(os.path.join(dst, "gross_region_density_boxplots_nc_tp.pdf"), bbox_inches = "tight")
plt.close()

#%%
#detailed structures
sois = ["Cerebral cortex", #telencephalon
        "Cerebral nuclei",
        "Thalamus",
        "Hypothalamus",
        "Ventral tegmental area", #diencephalon
        "Substantia nigra, compact part",
        "Substantia nigra, reticular part",
        "Red nucleus", #mesencephalon
        "Pontine gray", #metencephalon
        "Pontine reticular nucleus",
        "Inferior olivary complex"] #myelencephalon

#first calculate counts across entire region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    except:
        print(soi)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
    
counts_per_struct = np.array(counts_per_struct)

#voxels
vol = []
for soi in sois:
    progeny = []; counts = []; iids = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0])
    except:
        print(soi)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T)

# #boxplots of counts
# plt.figure(figsize = (5,5))
# df = pd.DataFrame(counts_per_struct.T)
# df.columns = sois
# g = sns.stripplot(data = df,  color = "dimgrey", orient = "h")
# sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
#             boxprops={'facecolor':'None'})

# plt.xlabel("Neurons")
# plt.ylabel("Region")
# plt.title("Neocortical timepoint")
# plt.savefig(os.path.join(dst, "detailed_region_counts_boxplots_nc_tp.pdf"), bbox_inches = "tight")
# plt.close()

#boxplots of density
plt.figure(figsize = (5,5))
df = pd.DataFrame(density)
df.columns = sois
g = sns.stripplot(data = df,  color = "#99C7E0", orient = "h")
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'})

plt.xlabel("Neurons/ mm$^3$")
plt.ylabel("Region")
plt.title("Neocortical timepoint")
plt.savefig(os.path.join(dst, "detailed_region_density_boxplots_nc_tp.pdf"), bbox_inches = "tight")
plt.close()