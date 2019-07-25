#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 16:11:23 2019

@author: wanglab
"""

from lightsheet.network_analysis import make_structure_objects
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl
from skimage.external import tifffile
import pandas as pd

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#custom
inj_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/injection_sites"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
cells_regions_pth = '/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/cell_count_by_region/nc_dataframe.p'

#making dictionary of injection sites
injections = {}

#MAKE SURE THEY ARE IN THIS ORDER
brains = ['20180409_jg46_bl6_lob6a_04', 
          '20180608_jg75', 
          '20170204_tp_bl6_cri_1750r_03',
          '20180608_jg72', 
          '20180416_jg56_bl6_lob8_04', 
          '20170116_tp_bl6_lob45_ml_11', 
          '20180417_jg60_bl6_cri_04', 
          '20180410_jg52_bl6_lob7_05', 
          '20170116_tp_bl6_lob7_1000r_10', 
          '20180409_jg44_bl6_lob6a_02', 
          '20180410_jg49_bl6_lob45_02', 
          '20180410_jg48_bl6_lob6a_01', 
          '20180612_jg80', 
          '20180608_jg71',
          '20170212_tp_bl6_crii_1000r_02', 
          '20170115_tp_bl6_lob6a_rpv_03', 
          '20170212_tp_bl6_crii_2000r_03', 
          '20180417_jg58_bl6_sim_02',
          '20170130_tp_bl6_sim_1750r_03', 
          '20170115_tp_bl6_lob6b_ml_04', 
          '20180410_jg50_bl6_lob6b_03', 
          '20170115_tp_bl6_lob6a_1000r_02', 
          '20170116_tp_bl6_lob45_500r_12', 
          '20180612_jg77', 
          '20180612_jg76', 
          '20180416_jg55_bl6_lob8_03', 
          '20170115_tp_bl6_lob6a_500r_01', 
          '20170130_tp_bl6_sim_rpv_01', 
          '20170204_tp_bl6_cri_1000r_02', 
          '20170212_tp_bl6_crii_250r_01', 
          '20180417_jg61_bl6_crii_05',
          '20170116_tp_bl6_lob7_ml_08', 
          '20180409_jg47_bl6_lob6a_05']

for pth in brains:
    print(pth)
    pth = os.path.join(inj_pth, pth+".tif.tif")
    injection = tifffile.imread(pth)
    print(injection.shape)
    injections[os.path.basename(pth)[:-8]] = injection #files have 2 .tif in the end
    
inj_raw = np.array([inj.astype(bool) for nm, inj in injections.items()])
    
atl_raw = tifffile.imread(atl_pth)
ann_raw = tifffile.imread(ann_pth)[:, 423:, :] #apparent cutoff
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

#%%

#pooled injections
ak_pool = np.asarray(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X",
                 "Simplex", "Crura", "PM, CP"])

#pooling injection regions
expr_all_as_frac_of_lob_pool = np.asarray([[brain[0]+brain[1]+brain[2]+brain[3], brain[4]+brain[5]+brain[6], 
                                 brain[7]+brain[8]+brain[9],brain[10], brain[11]+brain[12], 
                                 brain[13]+brain[14]] for brain in expr_all_as_frac_of_lob])

expr_all_as_frac_of_inj_pool = np.asarray([[brain[0]+brain[1]+brain[2]+brain[3], brain[4]+brain[5]+brain[6], 
                                 brain[7]+brain[8]+brain[9],brain[10], brain[11]+brain[12], 
                                 brain[13]+brain[14]] for brain in expr_all_as_frac_of_inj])
    
primary_pool = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj_pool])

#get n's after pooling
primary_lob_n = np.asarray([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])

#NORMALISATION
total_lob_sum = np.asarray([np.sum(expr_all_as_frac_of_lob_pool[i]) for i in range(len(expr_all_as_frac_of_lob))])    
normalised = np.asarray([expr_all_as_frac_of_lob_pool[i]/total_lob_sum[i] for i in range(len(expr_all_as_frac_of_lob))])

#%%
#making dictionary of cells by region
cells_regions = pckl.load(open(cells_regions_pth, "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

brains = list(cells_regions.keys())

#pooling regions
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#%%
#atlas res
scale_factor = 0.020 #mm/voxel

#GET ONLY STRIATUM small structures
sois = ["Caudoputamen","Nucleus accumbens", "Fundus of striatum", "Olfactory tubercle", 
        "Lateral septal nucleus", "Septofimbrial nucleus", "Septohippocampal nucleus",
        "Anterior amygdalar area", "Central amygdalar nucleus", "Intercalated amygdalar nucleus", "Medial amygdalar nucleus"]

#get densities for all the structures
df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", index_col = None)
df = df.drop(columns = ["Unnamed: 0"])
df = df.sort_values(by = ["name"])

#make new dict - for all brains
cells_pooled_regions = {} #for raw counts
volume_pooled_regions = {} #for density

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    print("\n\n"+brain)
    c_pooled_regions = {}
    d_pooled_regions = {}
    
    for soi in sois:
        print("\n\n"+soi)
        try:
            soi = [s for s in structures if s.name==soi][0]
            counts = [] #store counts in this list
            volume = [] #store volume in this list
            for k, v in cells_regions[brain].items():
                if k == soi.name:
                    counts.append(v)
            #add to volume list from LUT
            volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))
            progeny = [str(xx.name) for xx in soi.progeny]
            #now sum up progeny
            if len(progeny) > 0:
                for progen in progeny:
                    for k, v in cells_regions[brain].items():
                        if k == progen:
                            counts.append(v)
                            #add to volume list from LUT
                    volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0])
            c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
            d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)                    
            #add to volume list from LUT
            volume.append(df.loc[df.name == soi, "voxels_in_structure"].values[0])
            c_pooled_regions[soi] = np.sum(np.asarray(counts))
            d_pooled_regions[soi] = np.sum(np.asarray(volume))
                    
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions
    volume_pooled_regions[brain] = d_pooled_regions

    
#making the proper array per brain where regions are removed
cell_counts_per_brain = []

#initialise dummy var
i = []
for k,v in cells_pooled_regions.items():
    dct = cells_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    cell_counts_per_brain.append(i)
    #re-initialise for next
    i = []  
    
cell_counts_per_brain = np.asarray(cell_counts_per_brain)
#
#cell_counts_per_brain_pool = np.asarray([[brain[0]+brain[1]+brain[2]+brain[4], brain[3], brain[6], brain[5]+brain[7], 
#                                          brain[8]+brain[9], brain[10], brain[12], brain[11], 
#                                          brain[13]+brain[14], brain[15]+brain[16]] for brain in cell_counts_per_brain])

volume_per_brain = []

#initialise dummy var
i = []
for k,v in volume_pooled_regions.items():
    dct = volume_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    volume_per_brain.append(i)
    #re-initialise for next
    i = []  
#
#volume_per_brain_pool = np.asarray([[brain[0]+brain[1]+brain[2]+brain[4], brain[3], brain[6], brain[5]+brain[7], 
#                                          brain[8]+brain[9], brain[10], brain[12], brain[11], 
#                                          brain[13]+brain[14], brain[15]+brain[16]] for brain in volume_per_brain])
#    
volume_per_brain = np.asarray(volume_per_brain)*(scale_factor**3)
#volume_per_brain_pool = np.asarray(volume_per_brain_pool)*(scale_factor**3)

#calculate denisty
density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])
#density_per_brain_pool = np.asarray([xx/volume_per_brain_pool[i] for i, xx in enumerate(cell_counts_per_brain_pool)])

#%%
#get counts for all of neocortex
parent = ["Striatum"]    

#make new dict - for all brains
all_nc_counts = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    nc = {}
    
    for soi in parent:
        try:
            soi = [s for s in structures if s.name==soi][0]
            counts = [] #store counts in this list
            for k, v in cells_regions[brain].items():
                if k == soi.name:
                    counts.append(v)
            progeny = [str(xx.name) for xx in soi.progeny]
            #now sum up progeny
            if len(progeny) > 0:
                for progen in progeny:
                    for k, v in cells_regions[brain].items():
                        if k == progen:
                            counts.append(v)
            nc[soi.name] = np.sum(np.asarray(counts))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)
            nc[soi] = np.sum(np.asarray(counts))
                    
    #add to big dict
    all_nc_counts[brain] = nc
    
#making the proper array per brain where BRAIN NAMES are removed
nc_counts = []

#initialise dummy var
i = []
for k,v in all_nc_counts.items():
    dct = all_nc_counts[k]
    for j,l in dct.items():
        i.append(l)  
    nc_counts.append(i[0])
    #re-initialise for next
    i = []    
    
#make into % counts the proper way
cell_counts_per_brain_p = np.asarray([(brain/nc_counts[i]*100) for i, brain in enumerate(cell_counts_per_brain)])    
#%%    

#mean percent counts
mean_counts = np.asarray([np.mean(cell_counts_per_brain_p[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,4))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 10
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.3, aspect=10)
cb.set_label("% of striatal counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 1.5:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="xx-small")
dst = "/home/wanglab/Desktop"
plt.savefig(os.path.join(dst,"str_mean_percent_counts.pdf"), bbox_inches = "tight")

#%%

mean_counts = np.asarray([np.mean(density_per_brain[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 400
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,5)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.3, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 80:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="xx-small")
dst = "/home/wanglab/Desktop"
plt.savefig(os.path.join(dst,"str_mean_density.pdf"), bbox_inches = "tight")
#%%

#first, rearrange structures in ASCENDING order (will be plotted as descending) by density and counts
density_per_brain_descending_order = np.sort(density_per_brain)
order = np.argsort(np.mean(density_per_brain, axis = 0))
sois_descending_density = np.array(sois)[order]

cell_counts_per_brain_p_descending_order = np.sort(cell_counts_per_brain_p)
order = np.argsort(np.mean(cell_counts_per_brain_p, axis = 0))
sois_descending_pcounts = np.array(sois)[order]

#boxplots of densitys
plt.figure(figsize=(5,6))
plt.boxplot(density_per_brain_descending_order, vert = False, labels = sois_descending_density, sym = "", showcaps = False)
ngroup = len(density_per_brain_descending_order.T)
for i in range(ngroup):
    plt.scatter(density_per_brain_descending_order[:,i], y=np.ones(len(density_per_brain_descending_order[:,0]))*i+1, color = "k", s = 7)
plt.xlabel("Density ($cells/mm^3$)")
plt.ylabel("Striatal structures")
plt.savefig(os.path.join(dst,"str_density_boxplots.pdf"), bbox_inches = "tight")

#boxplots of percent counts
plt.figure(figsize=(5,6))
plt.boxplot(cell_counts_per_brain_p_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(cell_counts_per_brain_p_descending_order.T)
for i in range(ngroup):
    plt.scatter(cell_counts_per_brain_p_descending_order[:,i], 
                y=np.ones(len(cell_counts_per_brain_p_descending_order[:,0]))*i+1, color = "k", s = 10)
plt.xlabel("% of total striatal cells")
plt.ylabel("Striatal structures")
plt.savefig(os.path.join(dst,"str_pcounts_boxplots.pdf"), bbox_inches = "tight")

#%%

#GET ONLY STRIATUM small structures
sois = ["Globus pallidus, external segment", "Globus pallidus, internal segment", 
        "Substantia innominata", "Magnocellular nucleus", "Medial septal complex", 
        "Triangular nucleus of septum", "Bed nuclei of the stria terminalis", "Bed nucleus of the anterior commissure"]

#get densities for all the structures
df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", index_col = None)
df = df.drop(columns = ["Unnamed: 0"])
df = df.sort_values(by = ["name"])

#make new dict - for all brains
cells_pooled_regions = {} #for raw counts
volume_pooled_regions = {} #for density

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    print("\n\n"+brain)
    c_pooled_regions = {}
    d_pooled_regions = {}
    
    for soi in sois:
        print("\n\n"+soi)
        try:
            soi = [s for s in structures if s.name==soi][0]
            counts = [] #store counts in this list
            volume = [] #store volume in this list
            for k, v in cells_regions[brain].items():
                if k == soi.name:
                    counts.append(v)
            #add to volume list from LUT
            volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))
            progeny = [str(xx.name) for xx in soi.progeny]
            #now sum up progeny
            if len(progeny) > 0:
                for progen in progeny:
                    for k, v in cells_regions[brain].items():
                        if k == progen:
                            counts.append(v)
                            #add to volume list from LUT
                    volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0])
            c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
            d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)                    
            #add to volume list from LUT
            volume.append(df.loc[df.name == soi, "voxels_in_structure"].values[0])
            c_pooled_regions[soi] = np.sum(np.asarray(counts))
            d_pooled_regions[soi] = np.sum(np.asarray(volume))
                    
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions
    volume_pooled_regions[brain] = d_pooled_regions

    
#making the proper array per brain where regions are removed
cell_counts_per_brain = []

#initialise dummy var
i = []
for k,v in cells_pooled_regions.items():
    dct = cells_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    cell_counts_per_brain.append(i)
    #re-initialise for next
    i = []  
    
cell_counts_per_brain = np.asarray(cell_counts_per_brain)
#
#cell_counts_per_brain_pool = np.asarray([[brain[0]+brain[1]+brain[2]+brain[4], brain[3], brain[6], brain[5]+brain[7], 
#                                          brain[8]+brain[9], brain[10], brain[12], brain[11], 
#                                          brain[13]+brain[14], brain[15]+brain[16]] for brain in cell_counts_per_brain])

volume_per_brain = []

#initialise dummy var
i = []
for k,v in volume_pooled_regions.items():
    dct = volume_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    volume_per_brain.append(i)
    #re-initialise for next
    i = []  
#
#volume_per_brain_pool = np.asarray([[brain[0]+brain[1]+brain[2]+brain[4], brain[3], brain[6], brain[5]+brain[7], 
#                                          brain[8]+brain[9], brain[10], brain[12], brain[11], 
#                                          brain[13]+brain[14], brain[15]+brain[16]] for brain in volume_per_brain])
#    
volume_per_brain = np.asarray(volume_per_brain)*(scale_factor**3)
#volume_per_brain_pool = np.asarray(volume_per_brain_pool)*(scale_factor**3)

#calculate denisty
density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])
#density_per_brain_pool = np.asarray([xx/volume_per_brain_pool[i] for i, xx in enumerate(cell_counts_per_brain_pool)])

#%%
#get counts for all of neocortex
parent = ["Pallidum"]    

#make new dict - for all brains
all_nc_counts = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    nc = {}
    
    for soi in parent:
        try:
            soi = [s for s in structures if s.name==soi][0]
            counts = [] #store counts in this list
            for k, v in cells_regions[brain].items():
                if k == soi.name:
                    counts.append(v)
            progeny = [str(xx.name) for xx in soi.progeny]
            #now sum up progeny
            if len(progeny) > 0:
                for progen in progeny:
                    for k, v in cells_regions[brain].items():
                        if k == progen:
                            counts.append(v)
            nc[soi.name] = np.sum(np.asarray(counts))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)
            nc[soi] = np.sum(np.asarray(counts))
                    
    #add to big dict
    all_nc_counts[brain] = nc
    
#making the proper array per brain where BRAIN NAMES are removed
nc_counts = []

#initialise dummy var
i = []
for k,v in all_nc_counts.items():
    dct = all_nc_counts[k]
    for j,l in dct.items():
        i.append(l)  
    nc_counts.append(i[0])
    #re-initialise for next
    i = []    
    
#make into % counts the proper way
cell_counts_per_brain_p = np.asarray([(brain/nc_counts[i]*100) for i, brain in enumerate(cell_counts_per_brain)])    
#%%    

#mean percent counts
mean_counts = np.asarray([np.mean(cell_counts_per_brain_p[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(4,3))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 12
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,5)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.3, aspect=10)
cb.set_label("% of striatal counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 3:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="xx-small")
dst = "/home/wanglab/Desktop"
plt.savefig(os.path.join(dst,"pal_mean_percent_counts.pdf"), bbox_inches = "tight")

#%%

mean_counts = np.asarray([np.mean(density_per_brain[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(4,3))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 100
vmax = 800
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,8)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.3, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 200:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="xx-small")
dst = "/home/wanglab/Desktop"
plt.savefig(os.path.join(dst,"pal_mean_density.pdf"), bbox_inches = "tight")
#%%

#first, rearrange structures in ASCENDING order (will be plotted as descending) by density and counts
density_per_brain_descending_order = np.sort(density_per_brain)
order = np.argsort(np.mean(density_per_brain, axis = 0))
sois_descending_density = np.array(sois)[order]

cell_counts_per_brain_p_descending_order = np.sort(cell_counts_per_brain_p)
order = np.argsort(np.mean(cell_counts_per_brain_p, axis = 0))
sois_descending_pcounts = np.array(sois)[order]

#boxplots of densitys
plt.figure(figsize=(5,6))
plt.boxplot(density_per_brain_descending_order, vert = False, labels = sois_descending_density, sym = "", showcaps = False)
ngroup = len(density_per_brain_descending_order.T)
for i in range(ngroup):
    plt.scatter(density_per_brain_descending_order[:,i], y=np.ones(len(density_per_brain_descending_order[:,0]))*i+1, color = "k", s = 7)
plt.xlabel("Density ($cells/mm^3$)")
plt.ylabel("Pallidum structures")
plt.savefig(os.path.join(dst,"pal_density_boxplots.pdf"), bbox_inches = "tight")

#boxplots of percent counts
plt.figure(figsize=(5,6))
plt.boxplot(cell_counts_per_brain_p_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(cell_counts_per_brain_p_descending_order.T)
for i in range(ngroup):
    plt.scatter(cell_counts_per_brain_p_descending_order[:,i], 
                y=np.ones(len(cell_counts_per_brain_p_descending_order[:,0]))*i+1, color = "k", s = 10)
plt.xlabel("% of total pallidum cells")
plt.ylabel("Pallidum structures")
plt.savefig(os.path.join(dst,"pal_pcounts_boxplots.pdf"), bbox_inches = "tight")
