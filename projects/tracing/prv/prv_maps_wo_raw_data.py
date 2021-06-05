#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 19:04:53 2020

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, itertools, json, seaborn as sns, copy
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

#figure dest 
dst = "/home/wanglab/Desktop"

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
frac_of_inj_pool = data["frac_of_inj_pool"]
primary_pool = data["primary_pool"]
ak_pool = data["ak_pool"]
primary_lob_n = np.array([len(np.where(primary_pool == i)[0]) for i in range(max(primary_pool)+1)])

#change the lettering slightly 
ak_pool = np.array(["Lob. I-V", "Lob. VI, VII", "Lob. VIII-X",
       "Simplex", "Crus I", "Crus II", "PM, CP"])
#################################################SEPARATES INTO LAYER 5/6 CELLS, ONLY NEED TO RUN ONCE###########################################
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
with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for all of neocortex
sois = ["Frontal pole, cerebral cortex",
 "Infralimbic area",
 "Perirhinal area",
 "Gustatory areas",
 "Ectorhinal area",
 "Visceral area",
 "Prelimbic area",
 "Posterior parietal association areas",
 "Temporal association areas",
 "Orbital area",
 "Anterior cingulate area",
 "Auditory areas",
 "Agranular insular area",
 "Retrosplenial area",
 "Visual areas",
 "Somatomotor areas",
 "Somatosensory areas"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#then get only layers
layer56 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-7:] == "layer 5" or progen[-7:] == "Layer 5" or progen[-8:] == "layer 6a" or progen[-8:] == "layer 6b" or progen[-8:] == "Layer 6a" or progen[-8:] == "Layer 6b":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer56.append(np.array(counts).sum(axis = 0))
layer56 = np.array(layer56)        

#then get only layers
layer23 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-9:] == "layer 2/3" or progen[-9:] == "Layer 2/3":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer23.append(np.array(counts).sum(axis = 0))
layer23 = np.array(layer23)        

#calculate fraction of counts that come from layer 4,5,6, or 2/3
frac56 = np.sum(layer56, axis = 0)/np.sum(counts_per_struct, axis = 0)
mean_fracl56_per_struct = np.nanmean(frac56, axis = 0)

#layer p counts maps
pcounts = np.array([xx/sum(xx) for xx in counts_per_struct.T])*100

#get layer5/6 volumes
layer56_vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6a" or progen[-8:] == "layer 6b" or progen[-8:] == "Layer 6a" or progen[-8:] == "Layer 6b" or progen[-7:] == "layer 5" or progen[-7:] == "Layer 5":
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    layer56_vol.append(np.array(counts).sum(axis = 0))
layer56_vol = np.array(layer56_vol)        

density_l56 = np.array([xx/(layer56_vol[i]*(scale_factor**3)) for i, xx in enumerate(layer56)]).T
#%%
#clustering based on sam's MDS
cluster_num=[0]*5
cluster_num[0] = [0,6,14]
cluster_num[1] = [10,20,24]
cluster_num[2] = [3,8,9,11,12,19,21,23,23]
cluster_num[3] = [1,2,4,5,7,13,15,17]
cluster_num[4] = [16,18]
cluster_brains = [np.array(brains)[cluster_num[i]] for i in range(len(cluster_num))]
cluster_pcount = np.array([np.mean(pcounts[cluster_num[i]],axis=0) for i in range(len(cluster_num))]).T
cluster_inj = np.array([np.mean(frac_of_inj_pool[cluster_num[i]],axis=0) for i in range(len(cluster_num))]).T

#make % counts map 
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (1,4), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})
#inj fractions
ax = axes[0]
show = np.fliplr(cluster_inj)
cmap = copy.copy(plt.cm.RdPu)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.6
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")
ax.tick_params(length=6)

ax = axes[1]
show = cluster_pcount
# SET COLORMAP
vmin = 0
vmax = 10
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)
cb.set_label("mean % of\nneocortical neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
# aesthetics
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(sois, fontsize="x-small")
ax.set_xticks(np.arange(0, len(cluster_num))+.5)
ax.set_xticklabels(["A", "B", "C", "D", "E"], rotation=90)
plt.savefig(os.path.join(dst, "mds_mean_prv_pcounts_nc.svg"), bbox_inches = "tight")
plt.close()    
#%%
#ordered group heatmap
#cluster inds
cinds = list(itertools.chain.from_iterable(cluster_num))
#get pcounts per cluster
sort_pcounts = [pcounts[cluster_num[i]] for i in range(len(cluster_num))]
#get primary inj site per cluster
sort_pinj = [primary_pool[cluster_num[i]] for i in range(len(cluster_num))]
#sort primary inj site per cluster
sort_pcounts_ = [sort_pcounts[i][np.argsort(sort_pinj[i])] for i in range(len(sort_pcounts))]
#sort inj per cluster, and then by primary inj site per cluster
sort_inj = [frac_of_inj_pool[cluster_num[i]] for i in range(len(cluster_num))]
sort_inj_ = [sort_inj[i][np.argsort(sort_pinj[i])] for i in range(len(sort_inj))]
#load ml distances
from scipy.io import loadmat
mldist = loadmat("/jukebox/wang/zahra/tracing_projects/mapping_paper/count_data_final/prv_nc_mldist.mat")["inj"]
#flatten arr
mldist = np.array([xx[0] for xx in mldist])
#sort the same way by primary injection site
sort_ml = [mldist[cluster_num[i]] for i in range(len(cluster_num))]
sort_ml_ = [sort_ml[i][np.argsort(sort_pinj[i])] for i in range(len(sort_ml))]
#do for brains too..
sort_br = [np.array(brains)[cluster_num[i]] for i in range(len(cluster_num))]
sort_br_ = [sort_br[i][np.argsort(sort_pinj[i])] for i in range(len(sort_br))]
#flatten for maps
sort_pcounts__ = np.array(list(itertools.chain.from_iterable(sort_pcounts_)))
sort_inj__ = np.array(list(itertools.chain.from_iterable(sort_inj_)))
sort_ml__ = np.array(list(itertools.chain.from_iterable(sort_ml_)))
sort_br__ = np.array(list(itertools.chain.from_iterable(sort_br_)))
## display
fig, axes = plt.subplots(ncols = 1, nrows = 3, figsize = (6,12), 
                         sharex = False, gridspec_kw = {"wspace":0,"hspace":0,
                         "height_ratios": [1,0.15,8]})
#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj__).T
#SET COLORMAP HERE
cmap = copy.copy(plt.cm.RdPu)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, shrink=0.5, orientation="horizontal")#
cb.set_label("Injection % coverage\n of region", fontsize="small")
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.set_xticks([])
ax.tick_params(length=6)
#ml distances
ax = axes[1]
# ax.scatter(np.arange(len(sort_ml__)),sort_ml__,marker="|", LineWidth=2, color="k"); 
# ax.set_xlim([0,32])
# ax.axis("off")
# ax.set_yticks([70])
# ax.set_yticklabels(["Ml-distance"],fontsize="small")
pad = np.zeros((2,len(sort_ml__)))
pad[0] = sort_ml__
pad[1] = sort_ml__
show = np.absolute(pad) #absolute value bc the side of laterality doesn't matter
#colormap settings
cmap = copy.copy(plt.cm.Greens)
#colormap
pc = ax.pcolor(show, cmap=cmap)
cb = plt.colorbar(pc, shrink=0.5, orientation="horizontal")
cb.set_label("Medio-lateral distance (px)", fontsize="small")
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_xticks([])
ax.set_yticks([1])
ax.set_yticklabels(["Ml-distance"],fontsize="small")
ax.tick_params(length=6)
#pcounts
ax = axes[2]
show = np.flipud(np.fliplr(sort_pcounts__).T)
# SET COLORMAP
vmin = 0
vmax = 30
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, shrink=0.4, orientation="horizontal")
cb.set_label("% of neocortical neurons", fontsize="small")
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
# aesthetics
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(sois, fontsize="small")
ax.set_xticks([])
ax.set_xticklabels([])
ax.tick_params(length=6)
plt.savefig(os.path.join(dst, "mds_sorted_prv_pcounts_nc.svg"), bbox_inches="tight")

#%%
#hierarchical clustering
brains = np.array(brains)
df = pd.DataFrame(pcounts)
df.index = brains
df.columns = sois
#set cmap
maxpcount = 8
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0
vmax = maxpcount

h = sns.clustermap(df.T, cmap = cmap, row_cluster = False)
sns.despine(fig=None, ax=None, top=False, right=False, left=False, bottom=False, offset=None, trim=False)
plt.savefig(os.path.join(dst, "hierarchical_clustering_pcount_prv.svg"), bbox_inches="tight")
plt.close()

#order inj map by clusters
ind = h.dendrogram_col.reordered_ind
sort_inj = frac_of_inj_pool[ind]
sort_brains = brains[ind]
#make injection site heatmap only
fig, ax = plt.subplots(figsize = (5,2))
#inj fractions
show = np.fliplr(sort_inj).T
#colormap settings
cmap = copy.copy(plt.cm.Reds)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
lbls = np.asarray(sort_brains)
ax.set_xticklabels(lbls, rotation=90)
ax.tick_params(length=6)
plt.savefig(os.path.join(dst, "hierarchical_clustering_pcount_prv_inj.svg"), bbox_inches = "tight")
plt.close()   
#%%
#group into injection clusters based on drawn rectangles
cluster_num = [0, 3, 7, 10, 17, 21, 25]
sort_pcount = pcounts[ind]
cluster_brains = [sort_brains[cluster_num[i]:cluster_num[i+1]] for i in range(len(cluster_num)-1)]
cluster_pcount = np.array([np.mean(sort_pcount[cluster_num[i]:cluster_num[i+1]],axis=0) for i in range(len(cluster_num)-1)]).T
cluster_inj = np.array([np.mean(sort_inj[cluster_num[i]:cluster_num[i+1]],axis=0) for i in range(len(cluster_num)-1)]).T

#make % counts map 
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (1.3,4), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,5]})
#inj fractions
ax = axes[0]
show = np.fliplr(cluster_inj)
cmap = copy.copy(plt.cm.Reds)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.5
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")
ax.tick_params(length=6)

ax = axes[1]
show = cluster_pcount
# SET COLORMAP
vmin = 0
vmax = 11
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)
cb.set_label("% of neocortical neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
# aesthetics
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(sois, fontsize="x-small")
ax.set_xticks(np.arange(0, len(cluster_num)-1)+.5)
ax.set_xticklabels(np.arange(0, len(cluster_num)-1)+1)
plt.savefig(os.path.join(dst, "hclustering_mean_prv_pcounts_nc.svg"), bbox_inches = "tight")
plt.close()
#%%
#make injection site heatmap only
fig, ax = plt.subplots(figsize = (5,2))

sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
show = np.fliplr(sort_inj).T

cmap = plt.cm.RdPu 
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)#np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
lbls = np.asarray(sort_brains)
ax.set_xticklabels(np.array([ 1,  5, 10, 15, 20, 25]))
ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "prv_inj_nc.pdf"), bbox_inches = "tight")

#%%
#make % counts map like the h129 dataset (nc only for now)

## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (3.8,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})

#set colorbar features 
maxpcount = 30
whitetext = 3
yaxis = sois

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

# SET COLORMAP DETAILS HERE
cmap = plt.cm.RdPu 
cmap.set_under("white")
cmap.set_over(cmap(1.0))
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage \nof region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP

ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")

ax.tick_params(length=6)

ax = axes[1]
show = np.flipud(np.fliplr(sort_pcounts).T)

# SET COLORMAP
vmin = 0
vmax = maxpcount
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)
cb.set_label("% of neocortical neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# aesthetics
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(sois, fontsize="small")

ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "prv_pcounts_nc.pdf"), bbox_inches = "tight")

#%%

#make density map like the h129 dataset (nc only for now)

## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (3.8,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})

#set colorbar features 
maxdensity = 300
whitetext = 7
annotation_size = "small" #for the number annotations inside the heatmap

#sort inj fractions by primary lob
sort_density = [density_l56[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

cmap = plt.cm.RdPu 
cmap.set_under("white")
cmap.set_over(cmap(1.0))
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#, ticks=bounds, boundaries=bounds)
cb.set_label("Injection % coverage \nof region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.tick_params(length=6)

ax = axes[1]
show = np.flipud(np.fliplr(sort_density).T)

# SET COLORMAP
vmin = 0
vmax = maxdensity
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)
cb.set_label("Cells / mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="small")

ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "prv_density_nc.pdf"), bbox_inches = "tight")

#%%

#boxplots
#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
order = np.argsort(np.median(pcounts, axis = 0))[::-1]
#renaming for figure
sois_sort = np.array(sois)[order][:10]

#boxplots of percent counts
plt.figure(figsize = (5,4))
df = pd.DataFrame(pcounts)
df.columns = sois
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={"facecolor":"None"}, order = sois_sort)
plt.xlabel("% of neocortical neurons")
plt.ylabel("Region")

#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)

plt.tick_params(length=6)

plt.savefig(os.path.join(dst, "prv_nc_pcounts_boxplots.pdf"), bbox_inches = "tight")

#%%

#boxplots of density counts
order = np.argsort(np.median(density_l56, axis = 0))[::-1]

sois_sort = np.array(sois)[order][:10]

plt.figure(figsize = (5,4))
df = pd.DataFrame(density_l56)
df.columns = sois
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={"facecolor":"None"}, order = sois_sort)
plt.xlabel("Neurons / mm$^3$")
plt.ylabel("Region")

#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)
plt.tick_params(length=6)

plt.savefig(os.path.join(dst, "prv_nc_density_boxplots.pdf"), bbox_inches = "tight")

#%%

# SET COLORMAP
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#set min and max of colorbar
vmin = 0
vmax = 150

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(density_l56[np.where(primary_pool == idx)[0]], axis=0) 
    for idx in np.unique(primary_pool)])

fig, ax = plt.subplots(figsize=(2.8,6))

show = mean_counts.T 

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.3, aspect=10)
cb.set_label("Mean neurons / mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

#since we don't have a crus I primary injection site :(
ak_pool = np.array(["Lob. I-V", "Lob. VI, VII", "Lob. VIII-X", "Simplex",
       "Crus II", "PM, CP"])
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{} ({})".format(a, n) for a, n in zip(ak_pool, primary_lob_n[primary_lob_n > 0])], 
                    rotation = "vertical")

ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(sois)

plt.savefig(os.path.join(dst,"prv_nc_mean_density.pdf"), bbox_inches = "tight")