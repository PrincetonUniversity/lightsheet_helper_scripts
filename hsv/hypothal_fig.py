#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 17:07:32 2020

@author: zahra
"""

import matplotlib as mpl, json
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, pandas as pd, seaborn as sns

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

#from tom's script
brains = ["20170204_tp_bl6_cri_1000r_02", "20170115_tp_bl6_lob6a_rpv_03",
       "20170204_tp_bl6_cri_1750r_03", "20180417_jg58_bl6_sim_02",
       "20180410_jg51_bl6_lob6b_04", "20170204_tp_bl6_cri_250r_01",
       "20170130_tp_bl6_sim_rlat_05", "20180612_jg76",
       "20180417_jg57_bl6_sim_01", "20180409_jg45_bl6_lob6a_03",
       "20180417_jg60_bl6_cri_04", "20170130_tp_bl6_sim_1750r_03",
       "20170212_tp_bl6_crii_250r_01", "20170116_tp_bl6_lob45_500r_12",
       "20170212_tp_bl6_crii_1000r_02", "20170116_tp_bl6_lob45_ml_11",
       "20160916_tp_lob7_250r_05", "20180417_jg61_bl6_crii_05",
       "20170116_tp_bl6_lob7_1000r_10", "20180417_jg59_bl6_cri_03",
       "20170308_tp_bl6f_lob7_2x_02", "20170115_tp_bl6_lob6b_ml_04",
       "20180417_jg62_bl6_crii_06", "20180416_jg54_bl6_lob7_02",
       "20170115_tp_bl6_lob6a_1000r_02",
       "20170410_tp_bl6_lob6a_ml_repro_04", "20180410_jg50_bl6_lob6b_03",
       "20170207_db_bl6_crii_rpv_01", "20160916_tp_lob6_ml_04",
       "20161214_db_bl6_lob6a_lpvr_53d5hrs",
       "20161201_db_bl6_lob6b_500r_53d5hr"]

#figure dest 
dst = "/Users/zahra/Desktop"

###############################################################RUN AS IS#######################################################
#bucket path for data
src = "/Volumes/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_hypothalamus/dataframe_no_prog_at_each_level.p"
df_pth = "/Volumes/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ontology_file = "/Volumes/LightSheetTransfer/atlas/allen_atlas/allen.json"

cells_regions = pd.DataFrame(pckl.load(open(src, "rb"), encoding = "latin1"))
cells_regions["Structure"] = cells_regions.index

scale_factor = 0.020
try:
    ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])
except:
    ann_df = pd.read_excel(df_pth)

#%%
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

sois = ["Hypothalamus", 'Periventricular zone', 'Supraoptic nucleus',
       'Accessory supraoptic group',
       'Paraventricular hypothalamic nucleus',
       'Periventricular hypothalamic nucleus, anterior part',
       'Periventricular hypothalamic nucleus, intermediate part',
       'Arcuate hypothalamic nucleus', 'Periventricular region',
       'Anterodorsal preoptic nucleus', 'Anteroventral preoptic nucleus',
       'Anteroventral periventricular nucleus',
       'Dorsomedial nucleus of the hypothalamus',
       'Median preoptic nucleus', 'Medial preoptic area',
       'Vascular organ of the lamina terminalis',
       'Posterodorsal preoptic nucleus', 'Parastrial nucleus',
       'Periventricular hypothalamic nucleus, posterior part',
       'Periventricular hypothalamic nucleus, preoptic part',
       'Subparaventricular zone', 'Suprachiasmatic nucleus',
       'Subfornical organ', 'Ventrolateral preoptic nucleus',
       'Anterior hypothalamic nucleus', 'Mammillary body',
       'Lateral mammillary nucleus', 'Medial mammillary nucleus',
       'Supramammillary nucleus', 'Tuberomammillary nucleus',
       'Medial preoptic nucleus', 'Dorsal premammillary nucleus',
       'Ventral premammillary nucleus',
       'Ventromedial hypothalamic nucleus',
       'Posterior hypothalamic nucleus', 'Lateral hypothalamic area',
       'Lateral preoptic area', 'Preparasubthalamic nucleus',
       'Parasubthalamic nucleus', 'Retrochiasmatic area',
       'Subthalamic nucleus', 'Tuberal nucleus', 'Zona incerta',
       'Fields of Forel', 'Median eminence']

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

#drop striatum from rest
sois = sois[1:]

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

density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct[1:,:])]).T

#layer p counts maps
pcounts = np.array([xx[1:]/xx[0] for xx in counts_per_struct.T])*100

#sort pcounts and density by nuclei size
sois = np.array(sois)[np.argsort(vol)]
pcounts = pcounts.T[np.argsort(vol)].T
density = density.T[np.argsort(vol)].T

#%%
#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
order = np.argsort(np.median(density, axis = 0))[::-1]
sois_sort_density = np.array(sois)[order][:10]

order = np.argsort(np.median(pcounts, axis = 0))[::-1]
sois_sort_pcounts = np.array(sois)[order][:10]

#boxplots of densitys
plt.figure()
df = pd.DataFrame(density)
df.columns = sois 
g = sns.stripplot(data = df,  color = "#99C7E0", orient = "h", order = sois_sort_density)
sns.boxplot(data = df, orient = "h", showfliers=False,showcaps=False, boxprops={'facecolor':'None'}, order = sois_sort_density)
plt.xlabel("Neurons/ mm$^3$")
plt.savefig(os.path.join(dst, "hypothal_density_boxplots.pdf"), bbox_inches = "tight")

#boxplots of percent counts
plt.figure()
df = pd.DataFrame(pcounts)
df.columns = sois 
g = sns.stripplot(data = df,  color = "#99C7E0", orient = "h", order = sois_sort_pcounts)
sns.boxplot(data = df, orient = "h", showfliers=False,showcaps=False, boxprops={'facecolor':'None'}, order = sois_sort_pcounts)
plt.xlabel("% of total hypothalamus neurons")
plt.savefig(os.path.join(dst, "hypothal_pcounts_boxplots.pdf"), bbox_inches = "tight")


#boxplots of percent counts
plt.figure()
df = pd.DataFrame(pcounts)
df.columns = sois 
g = sns.stripplot(data = df,  color = "#99C7E0", orient = "h", order = sois_sort_pcounts)
sns.boxplot(data = df, orient = "h", showfliers=False,showcaps=False, boxprops={'facecolor':'None'}, order = sois_sort_pcounts)
plt.xlabel("% of total hypothalamus neurons")
plt.xlim([0, 40])
plt.savefig(os.path.join(dst, "hypothal_pcounts_boxplots_zoom.pdf"), bbox_inches = "tight")