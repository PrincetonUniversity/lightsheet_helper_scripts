#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 17:10:43 2020

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, itertools, json, seaborn as sns, copy
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl
import matplotlib

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
#sort inj fractions by primary lob
sort_pcounts = [pcounts[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_pcounts = np.array(list(itertools.chain.from_iterable(sort_pcounts)))
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))

#'Maximum projection percentage by lobule' - #TP EDITS
ncdata = np.copy(np.flipud(np.fliplr(sort_pcounts).T))

#threshold
threshold = 0
b = ncdata > threshold
ncdata = ncdata * b
#sns.heatmap(ncdata)

#now package into dataframe
ncdf = pd.DataFrame(np.flipud(ncdata), index=np.flipud(sois), columns=sort_brains).T

#add inj
idx_loc = {idx:loc for idx, loc in enumerate(ak_pool)}
nm_idx = {brain: idx for brain,idx in zip(brains, primary_pool)}
nm_loc = {nm:idx_loc[idx] for nm,idx in nm_idx.items()}
ncdf["brain"] = ncdf.index
ncdf["injection_loc"] = ncdf.apply(lambda xx: nm_loc[xx["brain"]],1)

#max proj
fig, ax = plt.subplots(figsize=(3,6))
ncdf = ncdf.groupby("injection_loc").max().drop(columns=["brain"])

#%%
#max pcount shuffle test
#make inj x nc region into an array
maxpcount = np.array(ncdf.T[ak_pool].T) #row = inj, column =  areas
X = maxpcount.T[::-1].T #note that this way, it is in the order of the sois (like the original pcount array)
#need to shuffle array by INJ 100 times
#compare to original maxpcount array? (t-test)
shuf = []
for itera in range(50):
    if itera%10 == 0: print(itera)        
        #rearrange inj sites
    s = X[np.random.choice(np.arange(len(X)), replace=False, size=len(X)),:]
    shuf.append(s)
shuf = np.array(shuf)
#quickly browse the shuffled data
plt.imshow(np.mean(shuf, axis=0).T, cmap = "Blues")
#test each of the % count observations (not just the max?) by the 100 max shuffles?
from scipy.stats import ttest_ind

pv_map = np.zeros(X.shape)
tstat_map = np.zeros(X.shape) #save both to figue out direction
for inj in range(X.shape[0]):
    print("\n"+ak_pool[inj])
    for region in range(X.shape[1]):
        b = shuf[:,inj,region] #somatosensory counts of max shuffle array
        a = pcounts[np.where(primary_pool==inj)][:,region] #somatosensory counts of all brains with a lob vi/vii injection
        s,pv = ttest_ind(a,b,equal_var=False) #inequal var does Welch’s t-test
        print("Region: %s\nP-value: %0.3f\nT-stat: %0.2f" % (sois[region],pv,s))
        pv_map[inj,region]=pv
        tstat_map[inj,region]=s

#run mutliple correction
from statsmodels.stats.multitest import multipletests
r,corrpv,alphasidak,alphabon = multipletests(np.ravel(pv_map),method="holm-sidak",alpha=1)
corrpv = corrpv.reshape(pv_map.shape)

#now display significance on top of max proj map
## display
fig,ax = plt.subplots(figsize=(3,6))
show = X.T 
# SET COLORMAP
vmin = 0
vmax = 30
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
annotation = False
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.3, aspect=10)
cb.set_label("Maximum % of neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
# signif
#only get positive significance           
sig = (corrpv.T < .05) & (tstat_map.T > 0)
for y,x in np.argwhere(sig):
    pass
    color="k"
    if show[y,x]>vmax-10:
        color="white"
    ax.text(x+0.5, y+0.4, "*", fontsize=12, horizontalalignment="center", verticalalignment="center",
            color=color, transform=ax.transData)
# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)
ax.set_xticklabels(ak_pool, rotation="vertical", fontsize=10)
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(sois, fontsize=10)
plt.title("*'s indicate signifance by Welch's t-test (p<0.05) and Holm-Sidak\n\
        multiple correction (p<0.1)\n", ha="center", va="center", fontsize="medium")
plt.savefig(os.path.join(fig_dst, "hsv_nc_maxpcount_shuffle.pdf"), bbox_inches = "tight")
plt.close()

#%%
#sort inj fractions by primary lob
sort_density = [density_l56[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))

#'Maximum projection DENSITY by lobule' - #TP EDITS
ncdata = np.copy(np.flipud(np.fliplr(sort_density).T))

#threshold
threshold = 0
b = ncdata > threshold
ncdata = ncdata * b
#sns.heatmap(ncdata)

#now package into dataframe
ncdf = pd.DataFrame(np.flipud(ncdata), index=np.flipud(sois), columns=sort_brains).T

#add inj
idx_loc = {idx:loc for idx, loc in enumerate(ak_pool)}
nm_idx = {brain: idx for brain,idx in zip(brains, primary_pool)}
nm_loc = {nm:idx_loc[idx] for nm,idx in nm_idx.items()}
ncdf["brain"] = ncdf.index
ncdf["injection_loc"] = ncdf.apply(lambda xx: nm_loc[xx["brain"]],1)
del ncdf["brain"]

#max proj
fig, ax = plt.subplots(figsize=(3,6))
ncdf = ncdf.groupby("injection_loc").max()

#%%
#max density shuffle test
#make inj x nc region into an array
maxd = np.array(ncdf) #row = inj, column = nc areas
X = maxd.T[::-1].T #note that this way, it is in the order of the sois (like the original pcount array)
#need to shuffle array by INJ some amount of times
#compare to original maxpcount array? (t-test)
shuf = []
for itera in range(20):
    # if itera%10 == 0: print(itera)        
    #rearrange inj sites
    s = X[np.random.choice(np.arange(len(X)), replace=False, size=len(X)),:]
    shuf.append(s)
shuf = np.array(shuf)
#quickly browse the shuffled data
plt.imshow(np.mean(shuf, axis=0).T, cmap = "Blues")
#test each of the % count observations (not just the max?) by the 100 max shuffles?
from scipy.stats import ttest_ind

pv_map = np.zeros(X.shape)
tstat_map = np.zeros(X.shape) #save both to figue out direction/tail
for inj in range(X.shape[0]):
    print("\n"+ak_pool[inj])
    for region in range(X.shape[1]):
        b = shuf[:,inj,region] #somatosensory counts of max shuffle array
        a = density_l56[np.where(primary_pool==inj)][:,region] #somatosensory counts of all brains with a lob vi/vii injection
        s,pv = ttest_ind(a,b,equal_var=False) #inequal var does Welch’s t-test
        print("No.of samples tested: %s\nRegion: %s\nP-value: %0.8f\nT-stat: %0.2f" % (len(b),sois[region],pv,s))
        pv_map[inj,region]=pv
        tstat_map[inj,region]=s

#run mutliple correction
#holm sidak
#cutoff=0.1
from statsmodels.stats.multitest import multipletests
r,corrpv,alphasidak,alphabon = multipletests(np.ravel(pv_map),
                                            method="holm-sidak",alpha=0.1)
# corrpv = corrpv.reshape(pv_map.shape)
corrpv=pv_map

#now display significance on top of max proj map
## display
fig,ax = plt.subplots(figsize=(3,6))
show = X.T 
# SET COLORMAP
vmin = 0
vmax = 500
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
annotation = False
#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, format="%d", shrink=0.3, aspect=10)
cb.set_label("Maximum density (Cells/$mm^3$)", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)
# signif
#only get positive significance           
sig = (corrpv.T < .05) & (tstat_map.T > 0)
for y,x in np.argwhere(sig):
    pass
    color="k"
    if show[y,x]>vmax-100:
        color="white"
    ax.text(x+0.5, y+0.4, "*", fontsize=12, horizontalalignment="center", verticalalignment="center",
            color=color, transform=ax.transData)
# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)
ax.set_xticklabels(ak_pool, rotation="vertical", fontsize=10)
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(sois, fontsize=10)
plt.title("*'s indicate signifance by Welch's t-test (p<0.05) and Holm-Sidak\n\
        multiple correction (p<0.1)\n", ha="center", va="center", fontsize="medium")
plt.savefig(os.path.join(fig_dst, "hsv_nc_maxdensity_shuffle.pdf"), bbox_inches = "tight")
plt.close()
