#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:53:28 2020

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, json, statsmodels.api as sm
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

#custom
src = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
cells_regions_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_contra_counts_33_brains_pma.csv"
dst = "/home/wanglab/Desktop/"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

data_pth = os.path.join(src, "nc_hsv_maps_contra_pma.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#set the appropritate variables
brains = data["brains"]
expr_all_as_frac_of_inj = data["expr_all_as_frac_of_inj"]
ak_pool = data["ak_pool"]
frac_of_lob = data["expr_all_as_frac_of_lob"]

brains = ['20180409_jg46_bl6_lob6a_04', '20180608_jg75',
       '20170204_tp_bl6_cri_1750r_03', '20180608_jg72',
       '20180416_jg56_bl6_lob8_04', '20170116_tp_bl6_lob45_ml_11',
       '20180417_jg60_bl6_cri_04', '20180410_jg52_bl6_lob7_05',
       '20170116_tp_bl6_lob7_1000r_10', '20180409_jg44_bl6_lob6a_02',
       '20180410_jg49_bl6_lob45_02', '20180410_jg48_bl6_lob6a_01',
       '20180612_jg80', '20180608_jg71', '20170212_tp_bl6_crii_1000r_02',
       '20170115_tp_bl6_lob6a_rpv_03', '20170212_tp_bl6_crii_2000r_03',
       '20180417_jg58_bl6_sim_02', '20170130_tp_bl6_sim_1750r_03',
       '20170115_tp_bl6_lob6b_ml_04', '20180410_jg50_bl6_lob6b_03',
       '20170115_tp_bl6_lob6a_1000r_02', '20170116_tp_bl6_lob45_500r_12',
       '20180612_jg77', '20180612_jg76', '20180416_jg55_bl6_lob8_03',
       '20170115_tp_bl6_lob6a_500r_01', '20170130_tp_bl6_sim_rpv_01',
       '20170204_tp_bl6_cri_1000r_02', '20170212_tp_bl6_crii_250r_01',
       '20180417_jg61_bl6_crii_05', '20170116_tp_bl6_lob7_ml_08',
       '20180409_jg47_bl6_lob6a_05']
#%%

cells_regions = pd.read_csv(cells_regions_pth)
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
sois = ['Somatosensory areas', 'Somatomotor areas', 'Visual areas',
       'Retrosplenial area', 'Agranular insular area', 'Auditory areas',
       'Anterior cingulate area', 'Orbital area',
       'Temporal association areas',
       'Posterior parietal association areas', 'Prelimbic area',
       'Visceral area', 'Ectorhinal area', 'Gustatory areas',
       'Perirhinal area', 'Infralimbic area',
       'Frontal pole, cerebral cortex']

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

layer_counts = np.array([layer1, layer23, layer4, layer5, layer6a, layer6b])

layers = np.array([np.mean(layer, axis = 1) for layer in layer_counts])

#%%

#make blue layer heatmap
fig, ax = plt.subplots(figsize = (2.7,6))

#inj fractions
show = np.flipud(layers.T)

cmap = plt.cm.Blues 
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 5
vmax = 180

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)#
cb.set_label("Mean neurons / mouse", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(sois))+.5)#np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(sois), fontsize="small")

ylbls = np.array([ "I", "II/III", "IV", "V", "VIa", "VIb"])
ax.set_xticks(np.arange(len(ylbls))+.5)
ax.set_xticklabels(ylbls)

plt.savefig(os.path.join(dst, "hsv_nc_layers.pdf"), bbox_inches = "tight")