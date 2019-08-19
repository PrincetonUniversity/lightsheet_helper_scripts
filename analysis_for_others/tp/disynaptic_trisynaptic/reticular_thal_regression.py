#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 14:59:54 2019

@author: wanglab
"""

import statsmodels.api as sm, os
import pickle, numpy as np, pandas as pd, matplotlib.pyplot as plt
from tools.analysis.network_analysis import make_structure_objects

nc_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/nc_dataframe_no_prog_at_each_level.p"
thal_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_thalamus/dataframe_no_prog_at_each_level.p"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

dst = "/home/wanglab/Desktop"
thal_data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data.p"
nc_data_pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
thal_data = pickle.load(open(thal_data_pth, "rb"), encoding = "latin1")
nc_data = pickle.load(open(nc_data_pth, "rb"), encoding = "latin1")


def get_counts_from_pickle(pth, sois, structures):
    cells_regions = pickle.load(open(pth, "rb"), encoding = "latin1")
   
    #get densities for all the structures
    df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", index_col = None)
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
    return brains, cell_counts_per_brain

thal_sois = ["Reticular nucleus of the thalamus"]

thal_brains, thal_counts_per_brain = get_counts_from_pickle(thal_pth, thal_sois, structures)
nc_brains, thalnc_counts_per_brain = get_counts_from_pickle(nc_pth, thal_sois, structures)

#GET ONLY NEOCORTICAL POOLS
nc_sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]
    
thal_brains, nct_counts_per_brain = get_counts_from_pickle(thal_pth, nc_sois, structures)
nc_brains, nc_counts_per_brain = get_counts_from_pickle(nc_pth, nc_sois, structures)

#sort vermis brains only
primary_pool = thal_data["primary_pool"]
thal_vermis = np.asarray([True if xx==0 or xx==1 else False for xx in primary_pool]) #mask

#thalamus
fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#linear regression
c = np.sum(nct_counts_per_brain, axis=1)[thal_vermis]
X = np.sort(c)
Y = np.squeeze(thal_counts_per_brain)[thal_vermis][np.argsort(c)]
Xlabels = np.array(thal_brains)[thal_vermis][np.argsort(c)]
results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
size = 40
ax.scatter(y = Y, x = X, s = size, color = "red")

ax.plot(mean_slope*range(4000)+mean_intercept, '--k')
    
lbls = Xlabels
#only show some labels 
#lbls = (np.zeros(Xlabels.shape)).astype('<U34')
#lbls[-5:] = Xlabels[-5:]
#lbls = ["" if xx=="0.0" else xx for xx in lbls]

for i, txt in enumerate(lbls):
    ax.annotate(txt, (X[i], Y[i]), fontsize = "x-small")
    
#ax.set_xlim([0, 100])
ax.set_xlabel("Total neocortical counts at thalamic timepoint")
ax.set_ylabel("Reticular nucleus counts")

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

plt.savefig(os.path.join(dst, "disynaptic.pdf"), bbox_inches = "tight")


#sort vermis brains only
primary_pool = nc_data["primary_pool"]
nc_vermis = np.asarray([True if xx==0 or xx==1 or xx==2 else False for xx in primary_pool]) #mask

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#neocortex
#linear regression
c = np.sum(nc_counts_per_brain, axis=1)[nc_vermis]
X = np.sort(c)
Y = np.squeeze(thalnc_counts_per_brain)[nc_vermis][np.argsort(c)]
Xlabels = np.array(nc_brains)[nc_vermis][np.argsort(c)]
results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
size = 40
ax.scatter(y = Y, x = X, s = size, color = "red")

ax.plot(mean_slope*range(30000)+mean_intercept, '--k')
    
#ax.set_xlim([0, 100])
ax.set_xlabel("Total neocortical counts at neocortical timepoint")
ax.set_ylabel("Reticular nucleus counts")

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

plt.savefig(os.path.join(dst, "disynaptic.pdf"), bbox_inches = "tight")
