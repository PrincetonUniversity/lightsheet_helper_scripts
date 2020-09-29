#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:45:59 2020

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, itertools, json, seaborn as sns, copy
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 4
mpl.rcParams["ytick.major.size"] = 4

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


#figure dest 
fig_dst = "/home/wanglab/Desktop"

###############################################################RUN AS IS#######################################################
#bucket path for data
src = "/jukebox/wang/zahra/tracing_projects/prv"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

cells_regions_pth = os.path.join(src, "for_tp/nc_contra_counts_25_brains_pma.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.020
try:
    ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])
except:
    ann_df = pd.read_excel(df_pth)

#imports
#path to pickle file
data_pth = os.path.join(src, "for_tp/prv_maps_contra_pma.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#set the appropritate variables
brains = data["brains"]
ak_pool = data["ak_pool"]
frac_of_inj_pool = data["frac_of_inj_pool"]
primary_pool = data["primary_pool"]
#get progeny of all large structures
with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for ss/sm areas
sois = ["Isocortex",
        "Primary somatosensory area, barrel field",
       "Primary somatosensory area, lower limb",
       "Primary somatosensory area, mouth",
       "Primary somatosensory area, nose",
       "Primary somatosensory area, trunk",
       "Primary somatosensory area, upper limb",
       "Supplemental somatosensory area",
       "Primary motor area", 
       "Secondary motor area"]

#first calculate counts across entire nc region
sssm_counts = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    sssm_counts.append(np.array(counts).sum(axis = 0))
sssm_counts = np.array(sssm_counts)

#get volumes
vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(sssm_counts)]).T #includes isocortex
#ss/sm p counts maps
pcounts = np.array([xx[1:]/xx[0] for ii,xx in enumerate(sssm_counts.T)])*100
#%%
#make % counts map like the h129 dataset (nc only for now)
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (8,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [3,5]})

#set colorbar features 
maxpcount = 20
yaxis = sois[1:] #for density by nc areas map

#sort inj fractions by primary lob
sort_pcounts = [pcounts[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_pcounts = np.array(list(itertools.chain.from_iterable(sort_pcounts)))
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

vmin = 0.05
vmax = 0.8
cmap = copy.copy(plt.cm.Reds)
cmap.set_over(cmap(1.0))
cmap.set_under("white")

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.5)
cb.set_label("Injection % coverage", fontsize="small")
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")

ax = axes[1]
show = np.fliplr(sort_pcounts).T

vmin = 0
vmax = maxpcount
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.5)
cb.set_label("% Neurons", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="small")

ax.set_xticks([])
ax.set_xticklabels([])

plt.savefig(os.path.join(fig_dst, "prv_pcounts_nc_sssm_areas.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(fig_dst, "prv_pcounts_nc_sssm_areas.jpg"), bbox_inches = "tight")
#%%

#density
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (8,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [3,5]})

#set colorbar features 
maxdensity = 300

#sort inj fractions by primary lob
sort_density = [density[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

vmin = 0.05
vmax = 0.8
cmap = copy.copy(plt.cm.Reds)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.5)
cb.set_label("Injection % coverage", fontsize="small")
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")

ax = axes[1]
show = np.fliplr(sort_density).T

vmin = 0
vmax = maxdensity
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
cmap.set_under("white")

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.5)
cb.set_label("Density (Cells/$mm^3$)", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="small")

ax.set_xticks([])
ax.set_xticklabels([])

plt.savefig(os.path.join(fig_dst, "prv_density_sssm_areas.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(fig_dst, "prv_density_sssm_areas.jpg"), bbox_inches = "tight")