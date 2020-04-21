#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:44:22 2019

@author: wanglab
"""

import matplotlib as mpl, itertools, pandas as pd, json
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#import data
pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/vtasnc_counts_contra.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")
cells_regions_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_contra_counts_23_brains_80um_ventric_erosion.csv"

#set dst 
dst = "/home/wanglab/Desktop"

frac_of_inj_pool = data["frac_of_inj_pool"]
brains = data["brainnames"]
ak_pool = data["ak_pool"]
primary_pool = data["primary_pool"]
short_nuclei = data["short_nuclei"]
primary_lob_n = data["primary_lob_n"]

#%%

#atlas res
scale_factor = 0.025 #mm/voxel

#GET ONLY VPM + bonus thal nuclei?, VTA, AND SNC COUNTS TO COMPARE
#rename structure column
df_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"

#get counts for all of neocortex
sois = ["Thalamus", "Midbrain",
        "Ventral tegmental area", #vta
        "Substantia nigra, reticular part", 
        "Substantia nigra, compact part",#snc
        "Reticular nucleus of the thalamus", #thal
        "Mediodorsal nucleus of thalamus",
        "Ventral posteromedial nucleus of the thalamus"
        ]

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.025
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
    if "msg" in list(dic.keys()): dic = dic["msg"][0]
    
    name = dic.get("name")
    children = dic.get("children")
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get("name")
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get("name")
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 

#get progeny of all large structures
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#then get volume
vol_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0]/2)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vol_per_struct.append(np.array(counts).sum(axis = 0))
vol_per_struct = np.array(vol_per_struct)        

density_per_struct = np.array([xx/(vol_per_struct[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)])[2:].T

short_nuclei = ["VTA", #vta
        "SNc", #snc
        "SNr", #snc
        "RTN", #thal
        "MD",
        "VPM"
        ]

#calculate pcounts by total cells in thalamus+midbrain?
pcounts = np.nan_to_num(np.array([(xx[2:]/(xx[0]+xx[1]))*100 for xx in counts_per_struct.T]))
#%%

#INJ FRACTION MAP

fig, ax = plt.subplots(figsize = (10,2))

sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
show = np.fliplr(sort_inj).T

vmin = 0.005
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over("darkred")
cmap.set_under("white")

#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("% coverage of lobule", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize="xx-small", ha="right")

#despline to make it look similar to paper figure
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.grid(False)

plt.savefig(os.path.join(dst, "vta_inj.pdf"), bbox_inches = "tight")

#%%
# SET COLORMAP
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))
annotations = False

#set min and max of colorbar
vmin = 0
vmax = 50

## mean density
fig = plt.figure(figsize=(5,2))
ax = fig.add_axes([.4,.1,.5,.8])

mean_d = np.asarray([np.mean(density_per_struct[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

show = mean_d.T #np.flip(mean_counts, axis = 1) # NOTE abs

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.8, aspect=12)

cb.set_label("Mean cells/$mm^3$", fontsize="small", labelpad=2)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
if annotations:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            pass
            if col > 0:
                ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="medium")
            if col > vmax-5:
                ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="medium")
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels([])#ak_pool, rotation="vertical", fontsize=8)
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)
ax.set_yticklabels([])#["{}".format(bi) for bi in short_nuclei], fontsize="medium")

#despline to make it look similar to paper figure
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.grid(False)
ax.tick_params(length=6)

plt.savefig(os.path.join(dst,"thalvtacomp_mean_density.pdf"), bbox_inches = "tight")
#%%

#set min and max of colorbar
vmin = 0
vmax = 4

## mean counts
fig = plt.figure(figsize=(5,2))
ax = fig.add_axes([.4,.1,.5,.8])

mean_counts = np.asarray([np.mean(pcounts[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

show = mean_counts.T 

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.8, aspect=12)

cb.set_label("Mean % of thalamic\nand midbrain neurons", fontsize="small", labelpad=8)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
if annotations:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            pass
            if col > 0:
                ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="medium")
            if col > vmax-5:
                ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="medium")
          
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels([])#ak_pool, rotation="vertical", fontsize=8)
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)
ax.set_yticklabels([])#["{}".format(bi) for bi in short_nuclei], fontsize="medium")

#despline to make it look similar to paper figure
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.grid(False)
ax.tick_params(length=6)

plt.savefig(os.path.join(dst,"thalvtacomp_mean_counts.pdf"), bbox_inches = "tight")