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

#NEOCORTICAL COHORT
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
dst = "/home/wanglab/Desktop"

#bucket path for data
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_hypothalamus/dataframe_no_prog_at_each_level.p"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

#need to segment and recalculate fractions
injsrc = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_injection_hypothalamus/hypothal_arrays"
inj_raw = np.array([tifffile.imread(os.path.join(injsrc,inj)).astype(bool) for inj in os.listdir(injsrc)])

atl_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif")[:, 423:, :] #cropping coz tom segmented this
ann_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")[:, 423:, :]
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)
     
#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lingula (I)": 912,
        "Lobule II": 976,
        "Lobule III": 984,
        "Lobule IV-V": 1091,
        "Lobule VIa": 936,
        "Lobule VIb": 1134,
        "Lobule VII": 944,
        "Lobule VIII": 951,
        "Lobule IX": 957, #uvula IX
        "Lobule X": 968, #nodulus X
        "Simplex lobule": 1007, #simplex
        "Crus 1": 1056, #crus 1
        "Crus 2": 1064, #crus 2
        "Paramedian lobule": 1025, #paramedian lob
        "Copula pyramidis": 1033 #copula pyramidis
        }
ak = np.asarray([k for k,v in iids.items()])

atlas_rois = {}
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    ann_blank = np.zeros_like(ann_raw)
    ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
    atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary

expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

#pooled injections
ak_pool = np.array(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X",
                 "Simplex", "Crus I", "Crus II", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
#get n's after pooling
primary_lob_n = np.array([len(np.where(primary_pool == i)[0]) for i in range(max(primary_pool)+1)])
#%%
cells_regions = pd.DataFrame(pckl.load(open(src, "rb"), encoding = "latin1"))
cells_regions["Structure"] = cells_regions.index

scale_factor = 0.020
try:
    ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])
except:
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
with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#only do top 20 nuclei by volume
nuclei = ["Hypothalamus",
 "Lateral hypothalamic area","Periventricular region",
 "Zona incerta","Mammillary body",
 "Anterior hypothalamic nucleus",
 "Periventricular zone","Lateral preoptic area",
 "Posterior hypothalamic nucleus",
 "Medial preoptic area","Tuberal nucleus",
 "Ventromedial hypothalamic nucleus",
 "Medial preoptic nucleus","Medial mammillary nucleus",
 "Dorsomedial nucleus of the hypothalamus",
 "Arcuate hypothalamic nucleus","Supramammillary nucleus",
 "Paraventricular hypothalamic nucleus",
 "Fields of Forel","Anteroventral periventricular nucleus",
 "Periventricular hypothalamic nucleus, intermediate part"]

#first calculate counts across entire region
counts_per_struct = []
for soi in nuclei:
    progeny = []; counts = []
    try:
        counts.append([cells_regions.loc[cells_regions.Structure==soi, brain].values[0] for brain in brains])
    except:
        counts.append([0]*len(brains))
        
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#get volumes
vol = []
for soi in nuclei:
    progeny = []; counts = []
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0])
    except:
        counts.append(0)
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

pcounts = np.nan_to_num(np.asarray([((brain[1:]/brain[0])*100) for brain in counts_per_struct.T]))    
density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T) #remove thalamus

#remove hypothalamus from density
sois = nuclei[1:]
density = density[:, 1:]
counts_per_struct = counts_per_struct[1:,:]
#%%
## display p counts
#set colorbar features 
maxpcount = 23
yaxis = np.flipud(sois)

#make % counts map like the h129 dataset
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (6.5,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,5]})

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
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.8)#
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
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.3)
cb.set_label("% of hypothalamic neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="small")
ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_pcounts_hypothal.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_pcounts_hypothal.jpg"), bbox_inches = "tight")
plt.close()
#%%
## display density

#set colorbar features 
maxd = 600
yaxis = np.flipud(sois)

#make density map like the h129 dataset
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (6.5,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,5]})


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
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.8)#
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
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.3)
cb.set_label("Neurons/ mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="small")
ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_density_hypothal.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(fig_dst, "nc_timepoint_hsv_density_hypothal.jpg"), bbox_inches = "tight")
plt.close()