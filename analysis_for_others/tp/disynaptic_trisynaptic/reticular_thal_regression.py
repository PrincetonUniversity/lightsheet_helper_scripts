#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 14:59:54 2019

@author: wanglab
"""

import statsmodels.api as sm, os
import pickle, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.ticker as ticker
from tools.analysis.network_analysis import make_structure_objects

nc_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/nc_dataframe_no_prog_at_each_level.p"
thal_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_thalamus/dataframe_no_prog_at_each_level.p"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

dst = "/home/wanglab/Desktop"
thal_data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data_v2.p"
nc_data_pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
thal_data = pickle.load(open(thal_data_pth, "rb"), encoding = "latin1")
nc_data = pickle.load(open(nc_data_pth, "rb"), encoding = "latin1")


def get_counts_from_pickle(pth, sois, structures, df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx",
                           scale = 0.020):
    cells_regions = pickle.load(open(pth, "rb"), encoding = "latin1")
   
    #get densities for all the structures
    df = pd.read_excel(df_pth, index_col = None)
    df = df.drop(columns = ["Unnamed: 0"])
    df = df.sort_values(by = ["name"])
    
    #make new dict - for all brains
    cells_pooled_regions = {} #for raw counts
    volume_pooled_regions = {} #for density
    
    brains = list(cells_regions.keys())
    
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
    i = []
    for k,v in volume_pooled_regions.items():
        dct = volume_pooled_regions[k]
        for j,l in dct.items():
            i.append(l)  
        volume_per_brain.append(i)
        #re-initialise for next
        i = []  
    volume_per_brain = np.asarray(volume_per_brain)*(scale**3)
    
    #calculate denisty
    density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])

    return brains, cell_counts_per_brain, density_per_brain

thal_sois = ["Reticular nucleus of the thalamus"]

thal_brains, thal_counts_per_brain, thal_density_per_brain = get_counts_from_pickle(thal_pth, thal_sois, structures)
nc_brains, thalnc_counts_per_brain, thalnc_density_per_brain = get_counts_from_pickle(nc_pth, thal_sois, structures)

#GET ONLY NEOCORTICAL POOLS
nc_sois = ["Isocortex"]
    
thal_brains, nct_counts_per_brain, nct_density_per_brain = get_counts_from_pickle(thal_pth, nc_sois, structures)
nc_brains, nc_counts_per_brain, nc_density_per_brain = get_counts_from_pickle(nc_pth, nc_sois, structures)

#mask to remove sandy brains
curated_brains = [False, True, True, False, False, False, True, False, True, False, True, True, False, False, 
                  False, False, True, False, True, False, True, False, True]

#sort vermis brains only
primary_pool = thal_data["primary_pool"][curated_brains]
thal_vermis = np.asarray([True if xx==0 or xx==1 or xx==3 else False for xx in primary_pool]) #mask
thal_density_per_brain = thal_density_per_brain[curated_brains]
nct_density_per_brain = nct_density_per_brain[curated_brains]
thal_brains = np.array(thal_brains)[curated_brains]
#sort nc vermis brains only
primary_pool = nc_data["primary_pool"]
nc_vermis = np.asarray([True if xx==0 or xx==1 or xx==2 or xx==4 else False for xx in primary_pool]) #mask for vermis and crus brains

#%%
vermis = True
#thalamus
fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

if vermis:
    #linear regression
    c = np.squeeze(nct_density_per_brain)[thal_vermis]
    X = np.sort(c)
    Y = np.squeeze(thal_density_per_brain)[thal_vermis][np.argsort(c)]
else:
    c = np.squeeze(nct_density_per_brain)
    X = np.sort(c)
    Y = np.squeeze(thal_density_per_brain)[np.argsort(c)]
    
Xlabels = np.array(thal_brains)[thal_vermis][np.argsort(c)]
results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
size = 70
ax.scatter(y = Y, x = X, s = size, color = "red")

#plot fit line
#ax.plot(mean_slope*range(50)+mean_intercept, '--k')
    
lbls = Xlabels
#only show some labels 
#lbls = (np.zeros(Xlabels.shape)).astype('<U34')
#lbls[-5:] = Xlabels[-5:]
#lbls = ["" if xx=="0.0" else xx for xx in lbls]

for i, txt in enumerate(lbls):
    ax.annotate(txt, (X[i], Y[i]), fontsize = "x-small")
    
ytick_spacing = 20; xtick_spacing = 2
ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
ax.set_xlim([0, 50])
ax.set_xlabel("Total neocortical density at thalamic timepoint")
ax.set_ylabel("Reticular nucleus density")

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

plt.savefig(os.path.join(dst, "thal_disynaptic.pdf"), bbox_inches = "tight")

#%%
#neocortex
#linear regression
if vermis:
    c = np.squeeze(nc_density_per_brain)[nc_vermis]
    X = np.sort(c)
    Y = np.squeeze(thalnc_density_per_brain)[nc_vermis][np.argsort(c)]
else:
    c = np.squeeze(nc_density_per_brain)
    X = np.sort(c)
    Y = np.squeeze(thalnc_density_per_brain)[np.argsort(c)]
    
Xlabels = np.array(nc_brains)[nc_vermis][np.argsort(c)]
results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#plot as scatter   
size = 40
ax.scatter(y = Y, x = X, s = size, color = "red")

ytick_spacing = 100; xtick_spacing = 20
ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
ax.set_xlim([0, 300])
ax.set_xlabel("Total neocortical density at neocortical timepoint")
ax.set_ylabel("Reticular nucleus density")

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

plt.savefig(os.path.join(dst, "nc_disynaptic.pdf"), bbox_inches = "tight")
