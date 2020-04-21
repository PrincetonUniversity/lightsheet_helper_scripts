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
from tools.imageprocessing.orientation import fix_orientation
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass
from scipy.stats import median_absolute_deviation as mad
from scipy.ndimage.measurements import center_of_mass

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#USING 60um edge erosion and 80 um ventricular erosion for NC, as edge seems to be the break bpoint. No real effect for ventricular so will keep the same
ann_pth = "/jukebox/wang/pisano/Python/atlas/stepwise_erosion/annotation_sagittal_atlas_20um_iso_60um_edge_erosion_80um_ventricular_erosion.tif"#"/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif"

#make structures
df_pth = "/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx"

dst = "/home/wanglab/Desktop"
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#%%

#GET A-P order for NC regions
import json
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

ann = tifffile.imread(ann_pth)
df = pd.read_excel(df_pth)
#get ids of all nc regions
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

primary_iid = []

for soi in sois:
    progeny = [soi]
    get_progeny(ontology_dict, soi, progeny)
    iids = [df.loc[df.name == progen, "id"].values[0] for progen in progeny]
    for iid in iids[1:]:
        z,y,x = np.where(ann == iid)
        ann[z,y,x] = iids[0]
    primary_iid.append(iids[0])
        

#reslice coronal
ann = fix_orientation(ann, ("2","0","1"))
ap_dim = []
for iid in primary_iid:
    ann_mask = ann.copy()
    ann_mask[ann_mask != iid] = 0
    ap, dv, ml = center_of_mass(ann_mask)
    ap_dim.append(ap)
    
ap_dim = np.array(ap_dim)
ap_dim_sort = np.sort(ap_dim)
sois = np.array(sois)[np.argsort(ap_dim)] #now use this ordering so it is truly ordered by A-P dimension
#%%
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

cells_regions_pth = "/home/wanglab/mounts/wang/zahra/h129_contra_vs_ipsi/data/thal_contra_counts_23_brains.csv"

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])

ann_df = "/home/wanglab/mounts/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
scale_factor = 0.025
ann_df = pd.read_excel(ann_df).drop(columns = ["Unnamed: 0"])

#make new dict - for all brains
cells_pooled_regions = {} #for raw counts
vol_pooled_regions = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    v_pooled_regions = {}
    
    for soi in sois:
        counts = []; vol = []
        soi = [s for s in structures if s.name==soi][0]
        
        counts.append(cells_regions.loc[cells_regions.Structure == soi.name, brain].values[0]) #store counts in this list
        vol.append(ann_df.loc[ann_df.name == soi.name, "voxels_in_structure"].values[0]) #store vols in this list
        #add to volume list from LUT
        progeny = [str(xx.name) for xx in soi.progeny]
        #now sum up progeny
        if len(progeny) > 0:
            for progen in progeny:
                counts.append(cells_regions.loc[cells_regions.Structure == progen, brain].values[0])
                if len(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values) > 0:
                    vol.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    
        c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
        v_pooled_regions[soi.name] = np.sum(np.asarray(vol))
    
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions
    vol_pooled_regions[brain] = v_pooled_regions
    
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

#making the proper array per brain where regions are removed
vol_per_brain = []

#initialise dummy var
i = []
for k,v in vol_pooled_regions.items():
    dct = vol_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    vol_per_brain.append(i)
    #re-initialise for next
    i = []
    
thal_cell_counts_per_brain = np.asarray(cell_counts_per_brain)
thal_vol_per_brain = np.array(vol_per_brain)/2
thal_density_per_brain = np.asarray([xx/(thal_vol_per_brain[i]*(scale_factor**3)) for i, xx in enumerate(thal_cell_counts_per_brain)])

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

cells_regions_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_contra_counts_33_brains_pma.csv"

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
ann_df = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
scale_factor = 0.020
ann_df = pd.read_excel(ann_df).drop(columns = ["Unnamed: 0"])

#make new dict - for all brains
cells_pooled_regions = {} #for raw counts
vol_pooled_regions = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    v_pooled_regions = {}
    
    for soi in sois:
        counts = []; vol = []
        soi = [s for s in structures if s.name==soi][0]
        
        counts.append(cells_regions.loc[cells_regions.Structure == soi.name, brain].values[0]) #store counts in this list
        vol.append(ann_df.loc[ann_df.name == soi.name, "voxels_in_structure"].values[0]) #store vols in this list
        #add to volume list from LUT
        progeny = [str(xx.name) for xx in soi.progeny]
        #now sum up progeny
        if len(progeny) > 0:
            for progen in progeny:
                counts.append(cells_regions.loc[cells_regions.Structure == progen, brain].values[0])
                vol.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    
        c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
        v_pooled_regions[soi.name] = np.sum(np.asarray(vol))
    
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions
    vol_pooled_regions[brain] = v_pooled_regions
    
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

#making the proper array per brain where regions are removed
vol_per_brain = []

#initialise dummy var
i = []
for k,v in vol_pooled_regions.items():
    dct = vol_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    vol_per_brain.append(i)
    #re-initialise for next
    i = []
    
nc_cell_counts_per_brain = np.asarray(cell_counts_per_brain)
nc_vol_per_brain = np.array(vol_per_brain)/2
nc_density_per_brain = np.asarray([xx/(nc_vol_per_brain[i]*(scale_factor**3)) for i, xx in enumerate(nc_cell_counts_per_brain)])
#%%
est_std = 0.6745 #normal mad to get estimated standard dev
mean_thal_counts = np.mean(thal_density_per_brain.sum(axis = 1), axis = 0)
mean_nc_counts = np.mean(nc_density_per_brain.sum(axis = 1), axis = 0)

ratio_mean = mean_thal_counts/mean_nc_counts

#calculate median also
median_thal_counts = np.median(thal_density_per_brain.sum(axis = 1), axis = 0)
median_nc_counts = np.median(nc_density_per_brain.sum(axis = 1), axis = 0)
ratio_median = median_thal_counts/median_nc_counts 

#st
ratio_std = np.std(thal_density_per_brain.sum(axis = 1), axis = 0)/np.std(nc_density_per_brain.sum(axis = 1), axis = 0)

#est std
mad_thal_counts = mad(thal_density_per_brain.sum(axis = 1), axis = 0)
mad_nc_counts = mad(nc_density_per_brain.sum(axis = 1), axis = 0)
ratio_mad = mad_thal_counts/mad_nc_counts
ratio_est_std = ratio_mad/est_std


#%%

import statsmodels.api as sm

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#rewrite nc labels
lbls = [df.loc[df.name == soi, "acronym"].values[0] for soi in sois]

#linear regression
Y = np.median(thal_density_per_brain, axis = 0)
X = range(len(lbls))
std = mad(thal_density_per_brain, axis = 0)

results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

ax.errorbar(X, Y, std, color='black', marker='o', linestyle='dashed', linewidth=2, markersize=8, label = "Slope=%0.2f\n$R$=%0.2f" % (mean_slope, np.sqrt(mean_r2)))

ax.set_ylabel("Neocortical density \nat thalamic timepoint \n(cells/$mm^3$)")
ax.set_xticks(np.arange(len(lbls))+.5)
ax.set_xticklabels(lbls, rotation=30, fontsize="small", ha="right")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")

ax.set_xlabel("$\Longleftrightarrow$ \n Anterior-Posterior")
plt.legend(prop={'size': 10}, bbox_to_anchor=(0.8,1), loc='upper left', ncol=1)

plt.savefig("/home/wanglab/Desktop/disynaptic_error_plot.pdf", bbox_inches = "tight")

#%%

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#rewrite nc labels
lbls = [df.loc[df.name == soi, "acronym"].values[0] for soi in sois]

#linear regression
Y = np.mean(thal_density_per_brain, axis = 0)
X = range(len(lbls))
std = np.std(thal_density_per_brain, axis = 0)

results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]
