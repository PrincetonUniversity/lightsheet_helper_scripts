#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 16:07:26 2020

@author: zahra
"""

import matplotlib as mpl, json
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, pandas as pd, seaborn as sns

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

#figure dest 
dst = "/home/wanglab/Desktop"

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
ak_pool = np.array(["Lobule I-V", "Lobule VI, VII", "Lobule VIII-X",
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

sois = ["Striatum", "Caudoputamen", "Nucleus accumbens",
       "Olfactory tubercle", "Lateral septal nucleus",
       "Medial amygdalar nucleus", "Central amygdalar nucleus",
       "Anterior amygdalar area", "Septofimbrial nucleus",
       "Fundus of striatum", "Intercalated amygdalar nucleus",
       "Septohippocampal nucleus"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    try:
        counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    except:
        counts.append([0]*len(brains))
        
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#get volumes
vol = []
for soi in sois:
    progeny = []; counts = []
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0]/2)
    except:
        counts.append(0)
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T

#layer p counts maps
pcounts = np.array([xx[1:]/xx[0] for xx in counts_per_struct.T])*100

#remove striatum from density
sois = sois[1:]
density = density[:, 1:]
counts_per_struct = counts_per_struct[1:,:]
#%%
## display p counts
#set colorbar features 
maxpcount = 30
yaxis = np.flipud(sois)

#make % counts map like the h129 dataset
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (6.5,3.5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [3,5]})

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

cmap = copy.copy(plt.cm.Reds)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.9)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.tick_params(length=6)

ax = axes[1]
show = np.fliplr(sort_pcounts).T

# SET COLORMAP
vmin = 0
vmax = maxpcount
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.7)
cb.set_label("% of striatum\nneurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="small")
ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_pcounts_striatum.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_pcounts_striatum.jpg"), bbox_inches = "tight")
plt.close()
#%%
## display density

#set colorbar features 
maxd = 800
yaxis = np.flipud(sois)

#make density map like the h129 dataset
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (6.5,3.5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [3,5]})

#sort inj fractions by primary lob
sort_density = [density[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

cmap = copy.copy(plt.cm.Reds)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.9)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.tick_params(length=6)

ax = axes[1]
show = np.fliplr(sort_density).T

# SET COLORMAP
vmin = 0
vmax = maxd
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.7)
cb.set_label("Neurons/ mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="small")
ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_density_striatum.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_density_striatum.jpg"), bbox_inches = "tight")
plt.close()

#%%
#mean percent counts
mean_counts = np.asarray([np.mean(pcounts[np.where(primary_pool == idx)[0]], 
                            axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,4))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T 

# SET COLORMAP
vmin = 0
vmax = 20
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
bounds = np.linspace(vmin,vmax,5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, shrink=0.5, aspect=10, format="%d")
cb.set_label("% of striatum neurons", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
        
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{} ({})".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation="vertical", fontsize="small")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="small")
plt.savefig(os.path.join(dst,"str_mean_pcounts.pdf"), bbox_inches = "tight")

#%%
mean_counts = np.asarray([np.mean(density[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,4))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T 

# SET COLORMAP
vmin = 0
vmax = 300
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, shrink=0.5, aspect=10, format="%d")
cb.set_label("Cells/ mm$^3$", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
        
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{} ({})".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation="vertical", fontsize="small")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="small")
plt.savefig(os.path.join(dst,"str_mean_density.pdf"), bbox_inches = "tight")

#%%
#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
order = np.argsort(np.median(density, axis = 0))[::-1]
sois_sort_density = np.array(sois)[order]

order = np.argsort(np.median(pcounts, axis = 0))[::-1]
sois_sort_pcounts = np.array(sois)[order]

#boxplots of densitys
plt.figure()
df = pd.DataFrame(density)
df.columns = sois 
g = sns.stripplot(data = df,  color = "steelblue", orient = "h", order = sois_sort_density)
sns.boxplot(data = df, orient = "h", showfliers=False,showcaps=False, boxprops={'facecolor':'None'}, order = sois_sort_density)
plt.xlabel("Cells/ mm$^3$")
plt.savefig(os.path.join(dst, "str_density_boxplots.pdf"), bbox_inches = "tight")

#boxplots of percent counts
plt.figure()
df = pd.DataFrame(pcounts)
df.columns = sois 
g = sns.stripplot(data = df,  color = "steelblue", orient = "h", order = sois_sort_pcounts)
sns.boxplot(data = df, orient = "h", showfliers=False,showcaps=False, boxprops={'facecolor':'None'}, order = sois_sort_pcounts)
plt.xlabel("% of total striatum neurons")
plt.savefig(os.path.join(dst, "str_pcounts_boxplots.pdf"), bbox_inches = "tight")