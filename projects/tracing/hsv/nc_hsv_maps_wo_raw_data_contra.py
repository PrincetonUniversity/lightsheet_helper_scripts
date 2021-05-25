#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 17:10:43 2020

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, itertools, json, seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

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
ak_pool = np.array(['Lob. I-V', 'Lob. VI, VII', 'Lob. VIII-X',
       'Simplex', 'Crus I', 'Crus II', 'PM, CP'])
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

#get counts for all of neocortex
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#get volumes
vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density_l56 = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T

#layer p counts maps
pcounts = np.array([xx/sum(xx) for xx in counts_per_struct.T])*100

#rename short sois
# sois = np.array(["IL", "PrL", "AC", "F Pole", "Orb", "Gust", "Insula", "Visc", "SM", "SS", "RS", "P Par", "VIS", 
#                     "Temp", "Aud", "EcR", "Pr"]) 

#sort pcounts and density by nuclei size
sois = np.array(sois)[np.argsort(vol)]
pcounts = pcounts.T[np.argsort(vol)].T
density_l56 = density_l56.T[np.argsort(vol)].T

#%%
#weighted sum of injection site
#idea is to make a boxplot/barplot of each region by injection site
#where the region is weighted by the injection fraction of that
#cerebellar region across brains

weighted_sum = [density_l56.T[i]*frac_of_inj_pool.T for i in range(len(sois))]
for i in range(len(sois)):
    plt.figure(figsize = (5,4))
    df = pd.DataFrame(weighted_sum[i].T)
    df.columns = ak_pool
    g = sns.stripplot(data = df,  color = "dimgrey", orient = "h")
    h = sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
                boxprops={"facecolor":"None"})
    plt.xlabel("%s\n\n cells / mm$^3$ x \nfraction of cerebellar region covered in injection" % sois[i])
    plt.ylabel("Cerebellar region")
    h.set_xlim([-0.5, 1000])
    h.set_xscale("symlog")
    plt.savefig(os.path.join(dst, "%s weighted_sum_density.pdf" % sois[i]), bbox_inches="tight")
    plt.close()

#pcounts
weighted_sum = [pcounts.T[i]*frac_of_inj_pool.T for i in range(len(sois))]
for i in range(len(sois)):
    plt.figure(figsize = (5,4))
    df = pd.DataFrame(weighted_sum[i].T)
    df.columns = ak_pool
    g = sns.stripplot(data = df,  color = "dimgrey", orient = "h")
    h = sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
                boxprops={"facecolor":"None"})
    plt.xlabel("%s\n\n cells / mm$^3$ x \nfraction of cerebellar region covered in injection" % sois[i])
    plt.ylabel("Cerebellar region")
    h.set_xlim([-0.5, 100])
    h.set_xscale("symlog")
    plt.savefig(os.path.join(dst, "%s weighted_sum_pcount.pdf" % sois[i]), bbox_inches="tight")
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

cmap = plt.cm.Oranges 
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
ax.set_xticklabels(np.array([ 1,  5, 10, 15, 20, 25, 30]))

ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "hsv_inj_orange_nc.pdf"), bbox_inches = "tight")

#%%
#set colorbar features 
maxpcount = 30
whitetext = 3
brain_lbl_size = "x-small"
yaxis = np.flipud(sois)


#make % counts map like the h129 dataset (nc only for now)
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (5,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})


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

#SET COLORMAP HERE
cmap = plt.cm.YlOrBr 
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
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")

#despline to make it look similar to paper figure
ax.tick_params(length=6)

ax = axes[1]
show = np.flipud(np.fliplr(sort_pcounts).T)

# SET COLORMAP
vmin = 0
vmax = maxpcount
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)
cb.set_label("% of neocortical neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="small")
ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

#despline to make it look similar to paper figure
ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "hsv_pcounts_nc_ylorbr_inj.pdf"), bbox_inches = "tight")


#%%

#make density map like the h129 dataset 
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (8,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})

#set colorbar features 
maxdensity = 400#300

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

cmap = plt.cm.YlOrBr
cmap.set_over(cmap(1.0))
vmin = 0.05
vmax = 0.8
cmap.set_under("white")

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#, ticks=bounds, boundaries=bounds)
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
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

#colorbar
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)
cb.set_label("Cells / mm$^3$", fontsize="x-small", labelpad=5)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)
# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="small")
ax.set_xticks(np.arange(0, len(sort_brains))+.5)
ax.set_xticklabels(sort_brains, rotation = "vertical")#np.arange(0, len(sort_brains), 5)+1)

ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "hsv_density_nc_ylorbr_inj.pdf"), bbox_inches = "tight")

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

plt.savefig(os.path.join(dst, "hsv_nc_pcounts_boxplots.pdf"), bbox_inches = "tight")

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
plt.xlabel("Cells / mm$^3$")
plt.ylabel("Region")

#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)

plt.tick_params(length=6)

plt.savefig(os.path.join(dst, "hsv_nc_density_boxplots.pdf"), bbox_inches = "tight")

#%%

# SET COLORMAP
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#set min and max of colorbar
vmin = 0
vmax = 400

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(density_l56[np.where(primary_pool == idx)[0]], axis=0) 
    for idx in np.unique(primary_pool)])

fig, ax = plt.subplots(figsize=(3,6))

show = mean_counts.T 

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.3, aspect=10)
cb.set_label("Mean neurons / mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
        
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{} ({})".format(a, n) for a, n in zip(ak_pool, primary_lob_n)], 
                    rotation = "vertical")

ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(sois)

plt.savefig(os.path.join(dst,"hsv_nc_mean_density.pdf"), bbox_inches = "tight")