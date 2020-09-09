#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 12:24:06 2020

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, itertools, json, seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl, copy 

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

#figure dest 
fig_dst = "/home/wanglab/Desktop"

###############################################################RUN AS IS#######################################################
#bucket path for data
src = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

cells_regions_pth = os.path.join(src, "nc_contra_counts_33_brains_pma.csv")

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
data_pth = os.path.join(src, "nc_hsv_maps_contra_pma.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#set the appropritate variables
brains = data["brains"]
expr_all_as_frac_of_inj = data["expr_all_as_frac_of_inj"]
ak_pool = data["ak_pool"]

#change the lettering slightly 
ak_pool = np.array(["Lob. I-V", "Lob. VI, VII", "Lob. VIII-X",
       "Simplex", "Crus I", "Crus II", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
primary_lob_n = np.array([len(np.where(primary_pool == i)[0]) for i in range(max(primary_pool)+1)])

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
with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for ss/sm areas
sois = ["Basolateral amygdalar nucleus", "Endopiriform nucleus",
        "Central amygdalar nucleus", "Cortical amygdalar area",
        "Lateral amygdalar nucleus", "Lateral amygdalar nucleus", 
        "Intercalated amygdalar nucleus", "Posterior amygdalar nucleus",
        "Anterior amygdalar area", "Somatomotor areas", "Somatosensory areas",
        "Auditory areas", "Visual areas", "Anterior cingulate area",
        "Retrosplenial area", "Temporal association areas", "Orbital area",
        "Visceral area", "Gustatory areas", "Posterior parietal association areas",
        "Ectorhinal area", "Prelimbic area", "Infralimbic area",
        "Perirhinal area", "Hippocampal region", "Retrohippocampal region",
        "Caudoputamen", "Lateral septal nucleus", "Pallidum",
        "Nucleus accumbens", "Fundus of striatum", "Septofimbrial nucleus",
        "Septohippocampal nucleus"]

#first calculate counts across entire nc region
str_counts = []
for soi in sois:
    # print(soi)
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add initial soi
    try:
        counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    except:
        counts.append([0]*len(brains))
    for progen in progeny:
        # print(progen)
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    str_counts.append(np.array(counts).sum(axis = 0))
str_counts = np.array(str_counts)

#get volumes
vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add initial soi
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0]/2)
    except:
        counts.append([0]*len(brains))
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vol.append(np.array(counts).sum(axis = 0))
    print(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(str_counts)]).T #includes isocortex
#p counts maps
#currently just summing across all the structures listed
pcounts = np.array([xx/sum(xx) for ii,xx in enumerate(str_counts.T)])*100

#only get counts for brains with lobule vi, crus 1, and crus 2 injections
#row =brain for pcounts array
pcounts_inj = pcounts[np.where((primary_pool==1) | (primary_pool==4) | (primary_pool==5))[0]]
density_inj = density[np.where((primary_pool==1) | (primary_pool==4) | (primary_pool==5))[0]]

#%%

#take mean of only lob vi,vii, crus 1, crus 2 injections
from matplotlib import colors 

cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))

#set min and max of colorbar
vmin = 0
vmax = 20

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(pcounts_inj[np.where(primary_pool[mask] == idx)[0]],axis=0) 
    for idx in [1,4,5]])

fig, ax = plt.subplots(figsize=(1.5,10))

show = np.flipud(mean_counts.T )

#colormap
pc = ax.pcolor(show, cmap=cmap, norm=colors.PowerNorm(gamma=0.4), vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Mean % neurons", fontsize="medium", labelpad=5)
cb.ax.tick_params(labelsize="medium")
cb.ax.set_visible(True)

#since we don't have a crus I primary injection site :(
ak_pool = np.array(["Lob. VI, VII", "Crus I", "Crus II"])
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{} ({})".format(a, n) for a, n in zip(ak_pool, [primary_lob_n[1],primary_lob_n[4],
                                                                     primary_lob_n[5]])], 
                    rotation = "vertical")

ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(np.array(sois)))

plt.savefig(os.path.join(fig_dst,"hsv_nc_mean_pcounts_jv.pdf"), bbox_inches = "tight")

#%%

#mean density
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))

#set min and max of colorbar
vmin = 0
vmax = 500

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(density_inj[np.where(primary_pool[mask] == idx)[0]],axis=0) 
    for idx in [1,4,5]])

fig, ax = plt.subplots(figsize=(1.5,10))

show = np.flipud(mean_counts.T )

#colormap
pc = ax.pcolor(show, cmap=cmap, norm=colors.PowerNorm(gamma=0.6), vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Mean neurons/$mm^3$", fontsize="medium", labelpad=5)
cb.ax.tick_params(labelsize="medium")
cb.ax.set_visible(True)

#since we don't have a crus I primary injection site :(
ak_pool = np.array(["Lob. VI, VII", "Crus I", "Crus II"])
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{} ({})".format(a, n) for a, n in zip(ak_pool, [primary_lob_n[1],primary_lob_n[4],
                                                                     primary_lob_n[5]])], 
                    rotation = "vertical")

ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(np.array(sois)))

plt.savefig(os.path.join(fig_dst,"hsv_nc_mean_density_jv.pdf"), bbox_inches = "tight")
#%%
#make % counts map like the h129 dataset (nc only for now)
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (8,8), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1,5]})

#set colorbar features 
maxpcount = 15
label_coordsy, label_coordsx  = -0.55,0.5 #for placement of vertical labels
brain_lbl_size = "x-small"
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

vmin = 0
vmax = 0.8
cmap = copy.copy(plt.cm.Reds)
cmap.set_over("darkred")
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.9, aspect=5)
cb.set_label("% of injection", fontsize="medium", labelpad=3)
cb.ax.tick_params(labelsize="medium")
cb.ax.set_visible(True)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="medium")

ax = axes[1]
show = np.fliplr(sort_pcounts).T

vmin = 0
vmax = maxpcount
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
#colormap
bounds = np.linspace(vmin,vmax,(int((vmax-vmin)/5)+1))

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("% Neurons", fontsize="medium", labelpad=3)
cb.ax.tick_params(labelsize="medium")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="medium")
ax.set_ylabel("Regions", fontsize="medium")
ax.yaxis.set_label_coords(label_coordsy, label_coordsx)

ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=brain_lbl_size, ha="right")

plt.savefig(os.path.join(fig_dst, "hsv_pcounts_nc_jess.pdf"), bbox_inches = "tight")
#%%
