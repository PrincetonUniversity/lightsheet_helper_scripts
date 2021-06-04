#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 18:48:40 2021

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl, json, itertools, statsmodels.api as sm
import tifffile, copy
from scipy.ndimage.measurements import center_of_mass
import matplotlib.colors
 
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"
fig_dst = "/home/wanglab/Desktop"
#pma atlas since it is just for injections
ann_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
#read annotation file
ann = tifffile.imread(ann_pth)
plt.imshow(ann[300])
z,y,x = ann.shape
#find 'medial' position (sagittal, so this will be the z value) 
medial = int(z/2)
#%%
#NEOCORTEX PRV
#brains
brains = ["20180205_jg_bl6f_prv_01", "20180205_jg_bl6f_prv_02", "20180205_jg_bl6f_prv_03", "20180205_jg_bl6f_prv_04", 
          "20180215_jg_bl6f_prv_05", "20180215_jg_bl6f_prv_06", "20180215_jg_bl6f_prv_08", "20180215_jg_bl6f_prv_09", 
           "20180305_jg_bl6f_prv_11", "20180305_jg_bl6f_prv_12", "20180305_jg_bl6f_prv_13","20180306_jg_bl6f_prv_14", 
           "20180305_jg_bl6f_prv_15", "20180312_jg_bl6f_prv_17", "20180326_jg_bl6f_prv_37",
           "20180313_jg_bl6f_prv_21", "20180313_jg_bl6f_prv_23", "20180313_jg_bl6f_prv_24", "20180313_jg_bl6f_prv_25",
           "20180322_jg_bl6f_prv_27", "20180322_jg_bl6f_prv_28", "20180323_jg_bl6f_prv_30", "20180326_jg_bl6f_prv_33", 
           "20180326_jg_bl6f_prv_34", "20180326_jg_bl6f_prv_35"]

src = "/jukebox/wang/zahra/tracing_projects/prv/prv_injection_sites"
imgs = [os.path.join(src, xx+".tif") for xx in brains]
#pool brain names and M-L dist into dict
prv_ml_dist = {}
#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    print(brain)
    inj_vol = tifffile.imread(img)
    z,y,x = inj_vol.shape
    z_c,y_c,x_c = center_of_mass(inj_vol)
    #take distance from center to arbitrary "midline" (aka half of z axis)
    dist = z_c-medial #positive = right, negative = left
    #save to dict 
    prv_ml_dist[brain[:-4]] = dist
#make heatmap of ML distances just to overlay on hclustering analysis
#(aka clustering by projection pattern)
#prv nc hclustering order
prv_brains = ['20180205_jg_bl6f_prv_01', '20180215_jg_bl6f_prv_08',
       '20180326_jg_bl6f_prv_37', '20180322_jg_bl6f_prv_28',
       '20180326_jg_bl6f_prv_35', '20180305_jg_bl6f_prv_11',
       '20180326_jg_bl6f_prv_34', '20180305_jg_bl6f_prv_12',
       '20180322_jg_bl6f_prv_27', '20180305_jg_bl6f_prv_13',
       '20180326_jg_bl6f_prv_33', '20180306_jg_bl6f_prv_14',
       '20180323_jg_bl6f_prv_30', '20180205_jg_bl6f_prv_04',
       '20180305_jg_bl6f_prv_15', '20180313_jg_bl6f_prv_23',
       '20180313_jg_bl6f_prv_25', '20180215_jg_bl6f_prv_06',
       '20180205_jg_bl6f_prv_03', '20180215_jg_bl6f_prv_09',
       '20180313_jg_bl6f_prv_21', '20180205_jg_bl6f_prv_02',
       '20180215_jg_bl6f_prv_05', '20180312_jg_bl6f_prv_17',
       '20180313_jg_bl6f_prv_24']
prv_sort_ml_dist = np.array([prv_ml_dist[brain] for brain in prv_brains]).T
#pad for plot
pad = np.zeros((2,len(prv_sort_ml_dist)))
pad[0] = prv_sort_ml_dist
pad[1] = prv_sort_ml_dist
#make ml distance heatmap only
fig, ax = plt.subplots(figsize = (6,0.5))
#inj fractions
show = np.absolute(pad) #absolute value bc the side of laterality doesn't matter
#colormap settings
cmap = copy.copy(plt.cm.Greens)
#colormap
pc = ax.pcolor(show, cmap=cmap)
cb = plt.colorbar(pc, shrink=5)
cb.set_label("Medio-lateral distance (px)", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_xticks([])
ax.set_yticks([])
ax.tick_params(length=6)
plt.savefig(os.path.join(fig_dst, "hclustering_ml_dist_nc_prv.svg"), bbox_inches = "tight")
plt.close()   
#%%
#NEOCORTEX HSV
#brains
brains = ["20170204_tp_bl6_cri_1750r_03",
       "20180416_jg56_bl6_lob8_04", "20170116_tp_bl6_lob45_ml_11",
       "20180410_jg52_bl6_lob7_05", "20170116_tp_bl6_lob7_1000r_10",
       "20180612_jg80", "20180608_jg71", "20170212_tp_bl6_crii_1000r_02",
       "20170115_tp_bl6_lob6a_rpv_03", "20170212_tp_bl6_crii_2000r_03",
       "20170130_tp_bl6_sim_1750r_03", "20170115_tp_bl6_lob6b_ml_04",
       "20170115_tp_bl6_lob6a_1000r_02", "20170116_tp_bl6_lob45_500r_12",
       "20180612_jg76", "20170115_tp_bl6_lob6a_500r_01",
       "20170130_tp_bl6_sim_rpv_01", "20170204_tp_bl6_cri_1000r_02",
       "20170212_tp_bl6_crii_250r_01", "20180409_jg46_bl6_lob6a_04",
       "20180608_jg75", "20180608_jg72", "20180417_jg60_bl6_cri_04",
       "20180409_jg44_bl6_lob6a_02", "20180410_jg49_bl6_lob45_02",
       "20180410_jg48_bl6_lob6a_01", "20180417_jg58_bl6_sim_02",
       "20180410_jg50_bl6_lob6b_03", "20180612_jg77",
       "20180416_jg55_bl6_lob8_03", "20180417_jg61_bl6_crii_05",
       "20170116_tp_bl6_lob7_ml_08", "20180409_jg47_bl6_lob6a_05"]
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_injection/nc_arrays/nc_arrays"
imgs = [os.path.join(src, xx+".tif.tif") for xx in brains]
#pool brain names and M-L dist into dict
nc_hsv_ml_dist = {}
#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    print(brain)
    inj_vol = tifffile.imread(img)
    z,y,x = inj_vol.shape
    z_c,y_c,x_c = center_of_mass(inj_vol)
    #take distance from center to arbitrary "midline" (aka half of z axis)
    dist = z_c-medial #positive = right, negative = left
    #save to dict 
    nc_hsv_ml_dist[brain[:-8]] = dist

#make heatmap of ML distances just to overlay on hclustering analysis
#(aka clustering by projection pattern)
#hsv nc hclustering order
nc_brains = ['20170115_tp_bl6_lob6a_1000r_02', '20170130_tp_bl6_sim_1750r_03',
       '20170116_tp_bl6_lob45_ml_11', '20170212_tp_bl6_crii_250r_01',
       '20180417_jg61_bl6_crii_05', '20170212_tp_bl6_crii_1000r_02',
       '20170116_tp_bl6_lob45_500r_12', '20170115_tp_bl6_lob6a_rpv_03',
       '20170115_tp_bl6_lob6b_ml_04', '20170204_tp_bl6_cri_1750r_03',
       '20180608_jg72', '20180612_jg77', '20180608_jg71',
       '20180409_jg46_bl6_lob6a_04', '20170130_tp_bl6_sim_rpv_01',
       '20180417_jg60_bl6_cri_04', '20180417_jg58_bl6_sim_02',
       '20180612_jg80', '20170212_tp_bl6_crii_2000r_03',
       '20180410_jg49_bl6_lob45_02', '20180410_jg52_bl6_lob7_05',
       '20170115_tp_bl6_lob6a_500r_01', '20180416_jg56_bl6_lob8_04',
       '20180410_jg48_bl6_lob6a_01', '20180612_jg76',
       '20170204_tp_bl6_cri_1000r_02', '20180608_jg75',
       '20170116_tp_bl6_lob7_1000r_10', '20180416_jg55_bl6_lob8_03',
       '20180409_jg44_bl6_lob6a_02', '20180410_jg50_bl6_lob6b_03',
       '20170116_tp_bl6_lob7_ml_08', '20180409_jg47_bl6_lob6a_05']

nc_hsv_sort_ml_dist = np.array([nc_hsv_ml_dist[brain] for brain in nc_brains]).T
#pad for plot
pad = np.zeros((2,len(nc_hsv_sort_ml_dist)))
pad[0] = nc_hsv_sort_ml_dist
pad[1] = nc_hsv_sort_ml_dist
#make ml distance heatmap only
fig, ax = plt.subplots(figsize = (6,0.5))
#inj fractions
show = np.absolute(pad)
#colormap settings
cmap = copy.copy(plt.cm.Greens)
#colormap
pc = ax.pcolor(show, cmap=cmap)
cb = plt.colorbar(pc, shrink=5)
cb.set_label("Medio-lateral distance (px)", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_xticks([])
ax.set_yticks([])
ax.tick_params(length=6)
plt.savefig(os.path.join(fig_dst, "hclustering_ml_dist_nc_hsv.svg"), bbox_inches = "tight")
plt.close()   
#%%
#THALAMUS HSV
#collect 
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20170410_tp_bl6_lob6a_ml_repro_01", "20160823_tp_bl6_cri_500r_02", "20180417_jg59_bl6_cri_03",
"20170207_db_bl6_crii_1300r_02", "20160622_db_bl6_unk_01", "20161205_tp_bl6_sim_750r_03",
"20180410_jg51_bl6_lob6b_04", "20170419_db_bl6_cri_rpv_53hr", "20170116_tp_bl6_lob6b_lpv_07",
"20170411_db_bl6_crii_mid_53hr", "20160822_tp_bl6_crii_1500r_06", "20160920_tp_bl6_lob7_500r_03",
"20170207_db_bl6_crii_rpv_01", "20161205_tp_bl6_sim_250r_02", "20161207_db_bl6_lob6a_500r_53hr",
"20170130_tp_bl6_sim_rlat_05", "20170115_tp_bl6_lob6b_500r_05", "20170419_db_bl6_cri_mid_53hr",
"20161207_db_bl6_lob6a_850r_53hr", "20160622_db_bl6_crii_52hr_01", "20161207_db_bl6_lob6a_50rml_53d5hr",
"20161205_tp_bl6_lob45_1000r_01", "20160801_db_l7_cri_01_mid_64hr"]    
#thalamus
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_injection_thalamus/201903_injection_thalamus/thal_arrays"
imgs = [os.path.join(src, xx+".tif.tif") for xx in brains]
#pool brain names and M-L dist into dict
thal_hsv_ml_dist = {}
#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    print(brain)
    inj_vol = tifffile.imread(img)
    z,y,x = inj_vol.shape
    z_c,y_c,x_c = center_of_mass(inj_vol)
    #take distance from center to arbitrary "midline" (aka half of z axis)
    dist = z_c-medial #positive = right, negative = left
    #save to dict 
    thal_hsv_ml_dist[brain[:-8]] = dist

#make heatmap of ML distances just to overlay on hclustering analysis
#(aka clustering by projection pattern)
#thalamus hclustering order
thal_brains = ['20160822_tp_bl6_crii_1500r_06', '20160622_db_bl6_crii_52hr_01',
       '20170115_tp_bl6_lob6b_500r_05', '20161205_tp_bl6_sim_250r_02',
       '20170411_db_bl6_crii_mid_53hr', '20160622_db_bl6_unk_01',
       '20161205_tp_bl6_sim_750r_03', '20160801_db_l7_cri_01_mid_64hr',
       '20160823_tp_bl6_cri_500r_02', '20170419_db_bl6_cri_rpv_53hr',
       '20170419_db_bl6_cri_mid_53hr',
       '20161207_db_bl6_lob6a_50rml_53d5hr',
       '20170410_tp_bl6_lob6a_ml_repro_01',
       '20161207_db_bl6_lob6a_850r_53hr', '20160920_tp_bl6_lob7_500r_03',
       '20170207_db_bl6_crii_rpv_01', '20161207_db_bl6_lob6a_500r_53hr',
       '20170116_tp_bl6_lob6b_lpv_07', '20170207_db_bl6_crii_1300r_02',
       '20170130_tp_bl6_sim_rlat_05', '20180417_jg59_bl6_cri_03',
       '20180410_jg51_bl6_lob6b_04', '20161205_tp_bl6_lob45_1000r_01']

thal_hsv_sort_ml_dist = np.array([thal_hsv_ml_dist[brain] for brain in thal_brains]).T
#pad for plot
pad = np.zeros((2,len(thal_hsv_sort_ml_dist)))
pad[0] = thal_hsv_sort_ml_dist
pad[1] = thal_hsv_sort_ml_dist
#make ml distance heatmap only
fig, ax = plt.subplots(figsize = (5,0.5))
#inj fractions
show = np.absolute(pad)
#colormap settings
cmap = copy.copy(plt.cm.Greens)
#colormap
pc = ax.pcolor(show, cmap=cmap)
cb = plt.colorbar(pc, shrink=5)
cb.set_label("Medio-lateral distance (px)", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_xticks([])
ax.set_yticks([])
ax.tick_params(length=6)
plt.savefig(os.path.join(fig_dst, "hclustering_ml_dist_thal_hsv.svg"), bbox_inches = "tight")
plt.close()   