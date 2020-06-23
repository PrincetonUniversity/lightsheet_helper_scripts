#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:22:57 2020

@author: wanglab
"""

import os, pandas as pd, numpy as np, pickle as pckl
import matplotlib.pyplot as plt, seaborn as sns, json, matplotlib as mpl
import itertools

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

#only get lobvi counts
lobvi = True
src = "/jukebox/wang/zahra/h129_contra_vs_ipsi"
df_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
cells_regions_pth_contra = os.path.join(src, "data/thal_contra_counts_23_brains_80um_ventric_erosion.csv")
cells_regions_pth_ipsi = os.path.join(src, "data/thal_ipsi_counts_23_brains_80um_ventric_erosion.csv")
dst = "/home/wanglab/Desktop"

#collect 
data_pth = os.path.join(src, "data/thal_hsv_maps_contra_allen.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

primary_pool = data["primary_pool"]
frac_of_inj_pool = data["frac_of_inj_pool"]
ak_pool = data["ak_pool"]
brains = np.array(data["brains"])
primary_lob_n = np.asarray([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])

#only get lobule VI brains
if lobvi: 
    brains = brains[primary_pool==1]
    frac_of_inj_pool = frac_of_inj_pool[primary_pool==1]

#get bilateral counts
#rearrange columns to match brain name
cells_regions_contra_w_structs = pd.read_csv(cells_regions_pth_contra)
cells_regions_contra = cells_regions_contra_w_structs[brains]
cells_regions_ipsi = pd.read_csv(cells_regions_pth_ipsi)[brains]
cells_regions = cells_regions_contra+cells_regions_ipsi
#add back structures column
cells_regions["Structure"] = cells_regions_contra_w_structs["Unnamed: 0"]
# cells_regions.to_csv(os.path.join(src, "data/thal_bilateral_counts_23_brains_80um_ventric_erosion.csv"))
scale_factor = 0.025
ann_df = pd.read_excel(df_pth)

def get_progeny(dic,parent_structure,progeny_list):
   
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

#sois = ["Cerebellar nuclei", "Vestibular nuclei", 
sois =  ["Pontine gray",
        "Tegmental reticular nucleus", "Lateral reticular nucleus",
        "Spinal nucleus of the trigeminal, caudal part",
        "Spinal nucleus of the trigeminal, interpolar part",
        "Spinal nucleus of the trigeminal, oral part",
        # "Principal sensory nucleus of the trigeminal",
        "External cuneate nucleus", "Gracile nucleus", "Cuneate nucleus"
        ]

#first calculate counts across entire region
counts_per_struct = []
for soi in sois:
    print(soi)
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

combined_sois = ["BPN/NRTP", "LRN", "Spinal trigeminal nuclei", 
                 "External cuneate nucleus", "Gracile nucleus", "Cuneate nucleus"]
combined_counts = np.array([np.array([xx[0]+xx[1], xx[2], xx[3]+xx[4]+xx[5],
                 xx[6],xx[7],xx[8]]) for xx in counts_per_struct.T]).T
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

combined_vol = np.array([vol[0]+vol[1], vol[2], vol[3]+vol[4]+vol[5],
                 vol[6],vol[7],vol[8]]).T
density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T) 
combined_density = np.nan_to_num(np.array([xx/(combined_vol[i]*(scale_factor**3)) for i, xx in enumerate(combined_counts)]).T) 

#%%
#display
#set colorbar features 
maxdensity = 300
#set true if need to sort structures in descending order of density/neurons
sort_descending = False

fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize = (5,6), sharex = False, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [5,5], "width_ratios": [15,1]})

#sort inj fractions by primary lob
if not lobvi:
    sort_density = [combined_density[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
    sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_brains = list(itertools.chain.from_iterable(sort_brains))
    sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))
else:
    sort_density = combined_density
    sort_brains = brains
    sort_inj = frac_of_inj_pool
    
if sort_descending:
    #now sort sois by # of neurons/density
    sort_sois = np.array(combined_sois)[np.argsort(np.median(sort_density,axis=0))]
    sort_density = sort_density.T[np.argsort(np.median(sort_density,axis=0))][::-1].T
    yaxis = sort_sois
else: 
    yaxis = np.flipud(combined_sois)
    
#inj fractions
ax = axes[0,0]
show = np.fliplr(sort_inj).T

cmap = plt.cm.Reds 
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="medium")
ax.tick_params(length=6)

ax = axes[1,0]
show = np.fliplr(sort_density).T

# SET COLORMAP
vmin = 0
vmax = maxdensity
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.6)
cb.set_label("Cells / mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="medium")
ax.set_xticks(np.arange(0, len(sort_brains))+.5)
ax.set_xticklabels(sort_brains, fontsize="x-small",rotation = "vertical")#np.arange(0, len(sort_brains), 5)+1)
ax.tick_params(length=6)

ax = axes[0,1]
ax.axis("off")
ax = axes[1,1]
show = np.flipud(np.array([np.mean(sort_density, axis=0)]).T)

# SET COLORMAP
vmin = 0
vmax = maxdensity
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)

# aesthetics
# yticks
ax.set_xticks(np.arange(1)+.5)
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels([])
ax.set_xticklabels(["Mean \ncells / mm$^3$"])#np.arange(0, len(sort_brains), 5)+1)
ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "lobvi_hsv_density_brainstem_dorsal_column_nuc_1_regions.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(dst, "lobvi_hsv_density_brainstem_dorsal_column_nuc_1_regions.jpg"), bbox_inches = "tight")

#%%
#display - just counts
#set colorbar features 
maxdensity = 400

#make density map like the h129 dataset
## display
fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize = (5,6), sharex = False, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [5,5], "width_ratios": [15,1]})

#sort inj fractions by primary lob
if not lobvi:
    sort_density = [combined_counts.T[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
    sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_brains = list(itertools.chain.from_iterable(sort_brains))
    sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))
else:
    sort_density = combined_counts.T
    sort_brains = brains
    sort_inj = frac_of_inj_pool
    
if sort_descending:
    #now sort sois by # of neurons/density
    sort_sois = np.array(combined_sois)[np.argsort(np.median(sort_density,axis=0))]
    sort_density = sort_density.T[np.argsort(np.median(sort_density,axis=0))][::-1].T
    yaxis = sort_sois
else: 
    yaxis = np.flipud(combined_sois)
    
#inj fractions
ax = axes[0,0]
show = np.fliplr(sort_inj).T

cmap = plt.cm.Reds 
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="medium")
ax.tick_params(length=6)

ax = axes[1,0]
show = np.fliplr(sort_density).T

# SET COLORMAP
vmin = 0
vmax = maxdensity
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.6)
cb.set_label("# Neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="medium")
ax.set_xticks(np.arange(0, len(sort_brains))+.5)
ax.set_xticklabels(sort_brains, fontsize="x-small",rotation = "vertical")#np.arange(0, len(sort_brains), 5)+1)
ax.tick_params(length=6)

ax = axes[0,1]
ax.axis("off")

ax = axes[1,1]
show = np.flipud(np.array([np.mean(sort_density, axis=0)]).T)

# SET COLORMAP
vmin = 0
vmax = maxdensity
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
# aesthetics
ax.set_xticks(np.arange(1)+.5)
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels([])
ax.set_xticklabels(["Mean \n# Neurons"])#np.arange(0, len(sort_brains), 5)+1)
ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "lobvi_hsv_counts_brainstem_dorsal_column_nuc_3_regions.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(dst, "lobvi_hsv_counts_brainstem_dorsal_column_nuc_3_regions.jpg"), bbox_inches = "tight")