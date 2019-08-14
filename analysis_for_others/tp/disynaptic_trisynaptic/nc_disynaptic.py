#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 21:00:17 2019

@author: wanglab
"""

%matplotlib inline
import numpy as np, pandas as pd, os, sys, shutil, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl
from tools.registration.register import elastix_command_line_call, jacobian_command_line_call, change_interpolation_order, transformix_command_line_call, transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.utils.io import listdirfull, makedir, load_memmap_arr, load_np, listall, load_kwargs
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#USING 60um edge erosion and 80 um ventricular erosion for NC, as edge seems to be the break bpoint. No real effect for ventricular so will keep the same
ann_pth = "/jukebox/wang/pisano/Python/atlas/stepwise_erosion/annotation_sagittal_atlas_20um_iso_60um_edge_erosion_80um_ventricular_erosion.tif"#"/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif"

#make structures
df_pth = "/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx"

dst = "/home/wanglab/Desktop"
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
structures_names = [xx.name for xx in structures]


#THALAMUS
brains = ["20170410_tp_bl6_lob6a_ml_repro_01",
         "20160823_tp_bl6_cri_500r_02",
         "20180417_jg59_bl6_cri_03",
         "20170207_db_bl6_crii_1300r_02",
         "20160622_db_bl6_unk_01",
         "20161205_tp_bl6_sim_750r_03",
         "20180410_jg51_bl6_lob6b_04",
         "20170419_db_bl6_cri_rpv_53hr",
         "20170116_tp_bl6_lob6b_lpv_07",
         "20170411_db_bl6_crii_mid_53hr",
         "20160822_tp_bl6_crii_1500r_06",
         "20160920_tp_bl6_lob7_500r_03",
         "20170207_db_bl6_crii_rpv_01",
         "20161205_tp_bl6_sim_250r_02",
         "20161207_db_bl6_lob6a_500r_53hr",
         "20170130_tp_bl6_sim_rlat_05",
         "20170115_tp_bl6_lob6b_500r_05",
         "20170419_db_bl6_cri_mid_53hr",
         "20161207_db_bl6_lob6a_850r_53hr",
         "20160622_db_bl6_crii_52hr_01",
         "20161207_db_bl6_lob6a_50rml_53d5hr",
         "20161205_tp_bl6_lob45_1000r_01",
         "20160801_db_l7_cri_01_mid_64hr"]

cells_regions_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_thalamus/dataframe_no_prog_at_each_level.p"

cells_regions = pckl.load(open(cells_regions_pth, "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

#get counts for all of neocortex
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#get densities for all the structures
df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", index_col = None)
df = df.drop(columns = ["Unnamed: 0"])
df = df.sort_values(by = ["name"])

#atlas res
scale_factor = 0.020 #mm/voxel

#make new dict - for all brains
cells_pooled_regions = {} #for raw counts
volume_pooled_regions = {} #for density

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    d_pooled_regions = {}
    
    for soi in sois:
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
    
volume_per_brain = np.asarray(volume_per_brain)*(scale_factor**3)

#calculate denisty
thal_density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])
mean_thal_density_per_brain = thal_density_per_brain.mean(axis = 0)
std_thal_density_per_brain = thal_density_per_brain.std(axis = 0)

#%%
## display
fig = plt.figure(figsize=(15, 5))
ax = fig.add_axes([.4,.1,.5,.8])

show = np.flipud(density_per_brain.T)

vmin = 0
vmax = 100
cmap = plt.cm.viridis
cmap.set_under("w")
cmap.set_over("gold")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.2, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 20:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# aesthetics
ax.set_xticks(np.arange(len(brains))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(brains)
ax.set_xticklabels(lbls, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(np.asarray(sois)), fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Thalamic 'disynaptic' timepoint", fontsize="x-small")
ax.yaxis.set_label_coords(-0.32,0.5)

plt.savefig(os.path.join(dst,"nc_density_at_thalamic_timepoint.pdf"), bbox_inches = "tight")

#%%
#NC
brains = ["20180409_jg46_bl6_lob6a_04",
 "20180608_jg75",
 "20170204_tp_bl6_cri_1750r_03",
 "20180608_jg72",
 "20180416_jg56_bl6_lob8_04",
 "20170116_tp_bl6_lob45_ml_11",
 "20180417_jg60_bl6_cri_04",
 "20180410_jg52_bl6_lob7_05",
 "20170116_tp_bl6_lob7_1000r_10",
 "20180409_jg44_bl6_lob6a_02",
 "20180410_jg49_bl6_lob45_02",
 "20180410_jg48_bl6_lob6a_01",
 "20180612_jg80",
 "20180608_jg71",
 "20170212_tp_bl6_crii_1000r_02",
 "20170115_tp_bl6_lob6a_rpv_03",
 "20170212_tp_bl6_crii_2000r_03",
 "20180417_jg58_bl6_sim_02",
 "20170130_tp_bl6_sim_1750r_03",
 "20170115_tp_bl6_lob6b_ml_04",
 "20180410_jg50_bl6_lob6b_03",
 "20170115_tp_bl6_lob6a_1000r_02",
 "20170116_tp_bl6_lob45_500r_12",
 "20180612_jg77",
 "20180612_jg76",
 "20180416_jg55_bl6_lob8_03",
 "20170115_tp_bl6_lob6a_500r_01",
 "20170130_tp_bl6_sim_rpv_01",
 "20170204_tp_bl6_cri_1000r_02",
 "20170212_tp_bl6_crii_250r_01",
 "20180417_jg61_bl6_crii_05",
 "20170116_tp_bl6_lob7_ml_08",
 "20180409_jg47_bl6_lob6a_05"]

cells_regions_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/nc_dataframe_no_prog_at_each_level.p"

cells_regions = pckl.load(open(cells_regions_pth, "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

#get densities for all the structures
df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", index_col = None)
df = df.drop(columns = ["Unnamed: 0"])
df = df.sort_values(by = ["name"])

#atlas res
scale_factor = 0.020 #mm/voxel

#make new dict - for all brains
cells_pooled_regions = {} #for raw counts
volume_pooled_regions = {} #for density

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    d_pooled_regions = {}
    
    for soi in sois:
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
    
volume_per_brain = np.asarray(volume_per_brain)*(scale_factor**3)

#calculate denisty
nc_density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])
mean_nc_density_per_brain = nc_density_per_brain.mean(axis = 0)
std_nc_density_per_brain = nc_density_per_brain.std(axis = 0)

#%%
## display
fig = plt.figure(figsize=(22, 5))
ax = fig.add_axes([.4,.1,.5,.8])

show = np.flipud(density_per_brain.T)

vmin = 0
vmax = 500
cmap = plt.cm.viridis
cmap.set_under("w")
cmap.set_over("gold")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.2, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 100:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(brains))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(brains)
ax.set_xticklabels(lbls, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(np.asarray(sois)), fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Neocortical 'trisynaptic' timepoint", fontsize="x-small")
ax.yaxis.set_label_coords(-0.22,0.5)

plt.savefig(os.path.join(dst,"nc_density_at_nc_timepoint.pdf"), bbox_inches = "tight")

#%%
ratio_mean_density = np.array(mean_thal_density_per_brain/mean_nc_density_per_brain)
ratio_std_density = np.array(std_thal_density_per_brain/std_nc_density_per_brain)

import pandas as pd
df = pd.DataFrame()
df["structures"] = sois
df["mean_thal_density"] = np.round(mean_thal_density_per_brain, decimals = 4)
df["mean_nc_density"] = np.round(mean_nc_density_per_brain, decimals = 4)
df["ratio_mean_density"] = np.round(ratio_mean_density, decimals = 4)
df["ratio_std_density"] = np.round(ratio_std_density, decimals = 4)

df.to_csv("/home/wanglab/Desktop/disynaptic.csv")

#%%

import statsmodels.api as sm

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#rewrite nc labels
lbls = ['Infralimbic',
 'Prelimbic',
 'Anterior cingulate',
 'Frontal pole',
 'Orbital',
 'Gustatory',
 'Agranular insula',
 'Visceral',
 'Somatosensory',
 'Somatomotor',
 'Retrosplenial',
 'Posterior parietal',
 'Visual',
 'Temporal',
 'Auditory',
 'Ectorhinal',
 'Perirhinal']
#linear regression
Y = np.sort(mean_nc_density_per_brain)[:-1]
X = mean_thal_density_per_brain[np.argsort(mean_nc_density_per_brain)][:-1]
strcs = np.asarray(lbls)[np.argsort(mean_nc_density_per_brain)][:-1]
results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[0]
mean_r2 = results.rsquared
mean_intercept = results.params[1]
#fit_thal = np.polyfit(range(thal_density_per_brain.shape[1]), thal_density_per_brain.mean(axis = 0), 1)
#fit_fn_thal = np.poly1d(fit_thal)
#linreg_stats_thal = linregress(range(thal_density_per_brain.shape[1]), thal_density_per_brain.mean(axis = 0))
#    
#plot as scatter
color_map = ["dimgray", "rosybrown", "darkred", "tomato", "chocolate", "orange", "gold", "olive", "darkseagreen", "springgreen", "teal",
             "darkturquoise", "steelblue", "navy", "indigo", "crimson", "deeppink"]

size = 30
for i in range(len(X)):
    ax.scatter(y = Y[i], x = X[i], s = size, label = strcs[i], c = color_map[i])
#
X = range(0, 30)
ax.plot(X, X*mean_slope + mean_intercept, "--r", label = "Slope={}\n$R^2$={}".format(round(mean_slope, 2), 
                   round(mean_r2, 2)))

##%%    
#plt.scatter(y = thal_density_per_brain[:, i], x = range(thal_density_per_brain.shape[0]), color = "red")
#    
ax.set_ylim([0, 500])
ax.set_xticks(np.arange(0, 30, 2))
#ax.set_xlim([0, 100])
ax.set_ylabel("Average neocortical densities at neocortical timepoint")
ax.set_xlabel("Average neocortical densities at thalamic timepoint")
plt.legend(prop={'size': 10}, bbox_to_anchor=(1,1), loc='upper left', ncol=1)
plt.savefig("/home/wanglab/Desktop/disynaptic.pdf", bbox_inches = "tight")
