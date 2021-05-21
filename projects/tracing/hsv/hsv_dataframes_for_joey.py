# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 15:13:12 2021

@author: SmartSPIM
"""

import os, pandas as pd, numpy as np, pickle as pckl
import matplotlib.pyplot as plt, seaborn as sns, json, matplotlib as mpl
import itertools

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

sois =  ["Cerebellar nuclei",
        "Inferior olivary complex"
        ]

#only get lobvi counts
lobvi = False
#change paths based on mounts
src =  r"W:\zahra\h129_contra_vs_ipsi"
df_pth = r"Z:\atlas\allen_atlas\allen_id_table_w_voxel_counts.xlsx" #"/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
cells_regions_pth_contra = os.path.join(src, r"data\thal_contra_counts_23_brains_80um_ventric_erosion.csv")
cells_regions_pth_ipsi = os.path.join(src, r"data\thal_ipsi_counts_23_brains_80um_ventric_erosion.csv")
dst = r"C:\Users\SmartSPIM\Desktop"#"/home/wanglab/Desktop"
#get progeny of all large structures
ontology_file = r"Z:\atlas\allen_atlas\allen.json" #"/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

#collect 
data_pth = os.path.join(src, r"data\thal_hsv_maps_contra_allen.p")
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

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)


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
#calculate density
density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T) 

#export to dataframe
#injection
df=pd.DataFrame(frac_of_inj_pool)
df.columns=ak_pool
df.index=brains
df.to_csv(r"Desktop\thal_inj_fractions.csv")
#COUNTS
df=pd.DataFrame(counts_per_struct.T)
df.index=brains
df.columns=["DCN_counts", "IO_counts"]
df["DCN_density"]=density[:,0]
df["IO_density"]=density[:,1]
df.to_csv(r"Desktop\thal_counts_density.csv")

#%%
#NEOCORTEX

#neocortex
src =  r"W:\zahra\h129_contra_vs_ipsi"
cells_regions_pth_contra = os.path.join(src, r"data\nc_contra_counts_33_brains_pma.csv")
cells_regions_pth_ipsi = os.path.join(src, r"data\nc_ipsi_counts_33_brains_pma.csv")

cells_regions_contra = pd.read_csv(cells_regions_pth_contra)
cells_regions_ipsi = pd.read_csv(cells_regions_pth_ipsi)
cells_regions_ipsi = cells_regions_ipsi[cells_regions_contra.columns]
cells_regions = cells_regions_ipsi + cells_regions_contra #bilateral

#rename structure column
cells_regions["Structure"] = cells_regions_contra["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])

scale_factor = 0.020
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

#collect
data_pth = os.path.join(src, r"data\nc_hsv_maps_contra_pma.p")
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
inj = frac_of_inj_pool #same as thal
primary = np.array([np.argmax(e) for e in frac_of_inj_pool])

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
#calculate density
density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T) 

#export to dataframe
#injection
df=pd.DataFrame(frac_of_inj_pool)
df.columns=ak_pool
df.index=brains
df.to_csv(r"Desktop\nc_inj_fractions.csv")
#COUNTS
df=pd.DataFrame(counts_per_struct.T)
df.index=brains
df.columns=["DCN_counts", "IO_counts"]
df["DCN_density"]=density[:,0]
df["IO_density"]=density[:,1]
df.to_csv(r"Desktop\nc_counts_density.csv")
#%%
#display
#set colorbar features 
maxdensity = 200
#set true if need to sort structures in descending order of density/neurons
sort_descending = False

fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize = (8,6), sharex = False, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [5,5], "width_ratios": [15,1]})

#divide density by maximum
norm_combined_density = np.array([xx/sum(xx) for xx in density])

#sort inj fractions by primary lob
if not lobvi:
    sort_density = [density[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
    sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_brains = list(itertools.chain.from_iterable(sort_brains))
    sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))
else:
    sort_density = density
    sort_brains = brains
    sort_inj = frac_of_inj_pool
    
if sort_descending:
    #now sort sois by # of neurons/density
    sort_sois = np.array(sois)[np.argsort(np.median(sort_density,axis=0))]
    sort_density = sort_density.T[np.argsort(np.median(sort_density,axis=0))][::-1].T
    yaxis = sort_sois
else: 
    yaxis = np.flipud(sois)
    
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
vmin = 0.1
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

plt.savefig(os.path.join(dst, "hsv_density_brainstem.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(dst, "hsv_density_brainstem.jpg"), bbox_inches = "tight")

#%%
#display - just counts
#set colorbar features 
maxdensity = 60

#make % counts array
pcounts = np.array([xx/sum(xx) for xx in counts.T])*100

#make density map like the h129 dataset
## display
fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize = (8,6), sharex = False, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [5,5], "width_ratios": [15,1]})

#sort inj fractions by primary lob
if not lobvi:
    sort_density = [pcounts[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
    sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
    sort_brains = list(itertools.chain.from_iterable(sort_brains))
    sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))
else:
    sort_density = pcounts
    sort_brains = brains
    sort_inj = frac_of_inj_pool
    
if sort_descending:
    #now sort sois by # of neurons/density
    sort_sois = np.array(sois)[np.argsort(np.median(sort_density,axis=0))]
    sort_density = sort_density.T[np.argsort(np.median(sort_density,axis=0))][::-1].T
    yaxis = sort_sois
else: 
    yaxis = np.flipud(sois)
    
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
vmin = 0.1
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

plt.savefig(os.path.join(dst, "hsv_pcounts_brainstem.pdf"), bbox_inches = "tight")
plt.savefig(os.path.join(dst, "hsv_pcounts_brainstem.jpg"), bbox_inches = "tight")