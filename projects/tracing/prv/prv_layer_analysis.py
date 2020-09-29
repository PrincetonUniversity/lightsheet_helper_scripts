#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:34:36 2020

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, json, statsmodels.api as sm, seaborn as sns, copy
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl


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

#custom
src = "/jukebox/wang/zahra/tracing_projects/prv"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
cells_regions_pth = os.path.join(src, "for_tp/nc_contra_counts_25_brains_pma.csv")
dst = "/home/wanglab/Desktop/"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

data_pth = os.path.join(src, "for_tp/prv_maps_contra_pma.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#set the appropritate variables
brains = data["brains"]
ak_pool = data["ak_pool"]
frac_of_inj_pool = data["frac_of_inj_pool"]

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.020
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

#get progeny of all large structures

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for all of neocortex

# sois =  ["Anterior cingulate area", "Orbital area",
#          "Prelimbic area","Infralimbic area",
#        "Frontal pole, cerebral cortex","Visual areas",
#        "Retrosplenial area", "Agranular insular area", "Auditory areas",
#        "Temporal association areas","Posterior parietal association areas",
#        "Visceral area", "Gustatory areas",
#        "Somatosensory areas", "Somatomotor areas","Ectorhinal area", "Perirhinal area","Entorhinal area"]

sois = ["Somatosensory areas", "Somatomotor areas", "Visual areas",
       "Retrosplenial area", "Agranular insular area", "Auditory areas",
       "Anterior cingulate area", "Orbital area",
       "Temporal association areas",
       "Posterior parietal association areas", "Prelimbic area",
       "Visceral area", "Ectorhinal area", "Gustatory areas",
       "Perirhinal area", "Infralimbic area",
       "Frontal pole, cerebral cortex"]

#get counts by layers
layer1 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-7:] == "layer 1" or progen[-7:] == "Layer 1":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer1.append(np.array(counts).sum(axis = 0))
layer1 = np.array(layer1)

layer23 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-9:] == "layer 2/3" or progen[-9:] == "Layer 2/3":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer23.append(np.array(counts).sum(axis = 0))
layer23 = np.array(layer23)

l4sois = ["Gustatory areas", "Visceral area", "Somatosensory areas", "Visual areas", "Temporal association areas",
            "Auditory areas"]
layer4 = []
for soi in sois:
    if soi not in l4sois:
        layer4.append([0]*len(brains))
    else:
        progeny = []; counts = []
        get_progeny(ontology_dict, soi, progeny)
        for progen in progeny:
            if progen[-7:] == "layer 4" or progen[-7:] == "Layer 4":
                counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
        layer4.append(np.array(counts).sum(axis = 0))
layer4 = np.array(layer4)

layer5 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-7:] == "layer 5" or progen[-7:] == "Layer 5":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer5.append(np.array(counts).sum(axis = 0))
layer5 = np.array(layer5)

layer6a = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6a" or progen[-8:] == "Layer 6a":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer6a.append(np.array(counts).sum(axis = 0))
layer6a = np.array(layer6a)

layer6b = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6b" or progen[-8:] == "Layer 6b":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer6b.append(np.array(counts).sum(axis = 0))
layer6b = np.array(layer6b)

#%%
#get vol
vollayer1 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-7:] == "layer 1" or progen[-7:] == "Layer 1":
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vollayer1.append(np.array(counts).sum(axis = 0))
vollayer1 = np.array(vollayer1)

vollayer23 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-9:] == "layer 2/3" or progen[-9:] == "Layer 2/3":
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vollayer23.append(np.array(counts).sum(axis = 0))
vollayer23 = np.array(vollayer23)

l4sois = ["Gustatory areas", "Visceral area", "Somatosensory areas", "Visual areas", "Temporal association areas",
            "Auditory areas"]
vollayer4 = []
for soi in sois:
    if soi not in l4sois:
        vollayer4.append(0)
    else:
        progeny = []; counts = []
        get_progeny(ontology_dict, soi, progeny)
        for progen in progeny:
            if progen[-7:] == "layer 4" or progen[-7:] == "Layer 4":
                counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
        vollayer4.append(np.array(counts).sum(axis = 0))
vollayer4 = np.array(vollayer4)

vollayer5 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-7:] == "layer 5" or progen[-7:] == "Layer 5":
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vollayer5.append(np.array(counts).sum(axis = 0))
vollayer5 = np.array(vollayer5)

vollayer6a = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6a" or progen[-8:] == "Layer 6a":
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vollayer6a.append(np.array(counts).sum(axis = 0))
vollayer6a = np.array(vollayer6a)

vollayer6b = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6b" or progen[-8:] == "Layer 6b":
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vollayer6b.append(np.array(counts).sum(axis = 0))
vollayer6b = np.array(vollayer6b)

layer_counts = np.array([layer1, layer23, layer4, layer5, layer6a, layer6b])
layers_counts_mean = np.array([np.mean(layer, axis = 1) for layer in layer_counts])

layer_vol = np.array([vollayer1, vollayer23, vollayer4, vollayer5, vollayer6a, vollayer6b])
layers_density = np.nan_to_num(np.array([np.array([xx/(layer_vol[l,i]*(scale_factor**3)) 
                                     for i, xx in enumerate(layer)]) 
                           for l,layer in enumerate(layer_counts)]).T)
#change inf numbers to 0
layers_density[layers_density < 1**100] = 0
layers_density_mean = layers_density.mean(axis = 0)
layers_density_mean[layers_density_mean == np.inf] = 0

#%%

#normalize p counts by region
n = np.nan_to_num(np.array([xx/sum(xx) for xx in layers_counts_mean.T]))*100

#make blue layer heatmap
fig, ax = plt.subplots(figsize = (2.7,6))
show = np.flipud(n)

cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 1
vmax = 40

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.4)#
cb.set_label("Mean % neurons\nper region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(sois), fontsize="small")

ylbls = np.array([ "I", "II/III", "IV", "V", "VIa", "VIb"])
ax.set_xticks(np.arange(len(ylbls))+.5)
ax.set_xticklabels(ylbls)

plt.savefig(os.path.join(dst, "prv_nc_layers_pcount_normalized.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(dst, "prv_nc_layers_pcount_normalized.jpg"), bbox_inches = "tight")
#%%

#make blue layer heatmap
fig, ax = plt.subplots(figsize = (2.7,6))
show = np.flipud(layers_density_mean)

cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 1
vmax = 300

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.4)#
cb.set_label("Mean neurons/ mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(sois), fontsize="small")

ylbls = np.array([ "I", "II/III", "IV", "V", "VIa", "VIb"])
ax.set_xticks(np.arange(len(ylbls))+.5)
ax.set_xticklabels(ylbls)

plt.savefig(os.path.join(dst, "prv_nc_layers_density.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(dst, "prv_nc_layers_density.jpg"), bbox_inches = "tight")

#%%
#now show only counts
#make blue layer heatmap
fig, ax = plt.subplots(figsize = (2.7,6))
show = np.flipud(layers_counts_mean.T)

cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 1
vmax = 200

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.4)#
cb.set_label("Mean cells/ mouse", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(sois), fontsize="small")

ylbls = np.array([ "I", "II/III", "IV", "V", "VIa", "VIb"])
ax.set_xticks(np.arange(len(ylbls))+.5)
ax.set_xticklabels(ylbls)

plt.savefig(os.path.join(dst, "hsv_nc_layers_counts.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(dst, "hsv_nc_layers_counts.jpg"), bbox_inches = "tight")

#%%

#make boxplots!!

#group into functional categories

func = np.array(["Sensory/motor", "Sensory/associative", "Frontal nonmotor", "Rhinal areas (memory)"])
markers = ["o", "v", "P", "D"]
flnms = ["sensmot", "sensassoc", "frontal", "rhinal"]
layers_counts_func = np.array([[[xx[0]+xx[1], xx[7]+xx[8]+xx[11]+xx[16], 
                                 xx[2]+xx[4]+xx[5]+xx[6]+xx[9]+xx[10]+xx[12]+xx[14]+xx[17],
                                 xx[3]+xx[15]+xx[13]] for xx in layer_counts[:, :, j]] for j in range(33)])
    
layers_pcounts_func = np.array([[xx/sum(xx) for xx in layers_counts_func[:, i, :]] for i in range(6)])*100

#rearrange for boxplots
t = layers_pcounts_func.transpose(1, 0, 2)
a = np.squeeze(np.array([[xx.ravel()] for xx in t])) #per area, per layer (e.g. area 1, layer 1, area 2, layer 1, etc...)

#%%
plt.figure(figsize = (5,10))

df = pd.DataFrame(a)
df.columns = ["{}, {}".format(f, l) for l in ylbls for f in func]
g = sns.stripplot(data = df,  color = "steelblue", orient = "h")
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={"facecolor":"None"})
plt.xlabel("% of neurons")
plt.ylabel("Region, Layer")

#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)

plt.savefig(os.path.join(dst, "hsv_nc_layers_pcount_boxplots.pdf"), bbox_inches = "tight")
 
#%%   
#get density boxplots
layer_vol_func = np.array([[xx[0]+xx[1], xx[7]+xx[8]+xx[11]+xx[16], 
                                 xx[2]+xx[4]+xx[5]+xx[6]+xx[9]+xx[10]+xx[12]+xx[14]+xx[17],
                                 xx[3]+xx[15]+xx[13]] for xx in layer_vol])   

layers_density_func = np.nan_to_num(np.array([np.array([xx/(layer_vol_func[l,i]*(scale_factor**3)) 
                                     for i, xx in enumerate(layer)]) 
                           for l,layer in enumerate(layers_counts_func.transpose(1, 2, 0))]))


#rearrange for boxplots
t = layers_density_func.transpose(2, 0, 1)
a = np.squeeze(np.array([[xx.ravel()] for xx in t])) #per area, per layer (e.g. area 1, layer 1, area 2, layer 1, etc...)

#%%
plt.figure(figsize = (5,10))

df = pd.DataFrame(a)
df.columns = ["{}, {}".format(f, l) for l in ylbls for f in func]
g = sns.stripplot(data = df,  color = "steelblue", orient = "h")
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={"facecolor":"None"})
plt.xlabel("Neurons/ mm$^3$")
plt.ylabel("Region, Layer")

g.set(xscale="log")

#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)

plt.savefig(os.path.join(dst, "hsv_nc_layers_density_boxplots.pdf"), bbox_inches = "tight")