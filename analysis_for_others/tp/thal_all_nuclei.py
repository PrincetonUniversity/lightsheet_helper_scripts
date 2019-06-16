#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 13:38:57 2019

@author: wanglab
"""

import sys
sys.path.append("/jukebox/wang/zahra/python/lightsheet_helper_scripts")
from lightsheet.network_analysis import make_structure_objects
from scipy.stats import spearmanr
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, pandas as pd
from skimage.external import tifffile
import seaborn as sns

#custom
inj_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
cells_regions_pth = '/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/cell_count_by_region/dataframe.p'

#making dictionary of injection sites
injections = {}

#MAKE SURE THEY ARE IN THIS ORDER
brains = ['20170410_tp_bl6_lob6a_ml_repro_01',
         '20160823_tp_bl6_cri_500r_02',
         '20180417_jg59_bl6_cri_03',
         '20170207_db_bl6_crii_1300r_02',
         '20160622_db_bl6_unk_01',
         '20161205_tp_bl6_sim_750r_03',
         '20180410_jg51_bl6_lob6b_04',
         '20170419_db_bl6_cri_rpv_53hr',
         '20170116_tp_bl6_lob6b_lpv_07',
         '20170411_db_bl6_crii_mid_53hr',
         '20160822_tp_bl6_crii_1500r_06',
         '20160920_tp_bl6_lob7_500r_03',
         '20170207_db_bl6_crii_rpv_01',
         '20161205_tp_bl6_sim_250r_02',
         '20161207_db_bl6_lob6a_500r_53hr',
         '20170130_tp_bl6_sim_rlat_05',
         '20170115_tp_bl6_lob6b_500r_05',
         '20170419_db_bl6_cri_mid_53hr',
         '20161207_db_bl6_lob6a_850r_53hr',
         '20160622_db_bl6_crii_52hr_01',
         '20161207_db_bl6_lob6a_50rml_53d5hr',
         '20161205_tp_bl6_lob45_1000r_01',
         '20160801_db_l7_cri_01_mid_64hr']

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
iids = {"Lobule III": 984,
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

#just cerebellar region names
ak = np.asarray([k for k,v in iids.items()])

atlas_rois = {}
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    ann_blank = np.zeros_like(ann_raw)
    ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
    atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary
    
#get fractions
expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

print(expr_all_as_frac_of_lob[15])

#%%
#making dictionary of cells by region
cells_regions = pckl.load(open(cells_regions_pth, "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      
#make sure brains are the same order
assert brains == list(cells_regions.keys())

#pooling regions
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#%%
#GET ONLY THALAMUS POOLS
nuclei = ["Parafascicular nucleus",
        "Posterior complex of the thalamus",        
        "Posterior triangular thalamic nucleus",
        "Lateral posterior nucleus of the thalamus",
        "Lateral habenula",
        "Lateral dorsal nucleus of thalamus",
        "Central lateral nucleus of the thalamus",        
        "Paraventricular nucleus of the thalamus",
        "Nucleus of reuniens",        
        "Mediodorsal nucleus of thalamus",
        "Ventral part of the lateral geniculate complex",
        "Ventral posterolateral nucleus of the thalamus",
        "Ventral posteromedial nucleus of the thalamus",         
        "Submedial nucleus of the thalamus",
        "Reticular nucleus of the thalamus",
        "Ventral medial nucleus of the thalamus",
        "Anteroventral nucleus of thalamus",
        "Ventral anterior-lateral complex of the thalamus"
        ]

#make new dict - for all brains
cells_pooled_regions = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    pooled_regions = {}
    
    for soi in nuclei:
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
            pooled_regions[soi.name] = np.sum(np.asarray(counts))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)
            pooled_regions[soi] = np.sum(np.asarray(counts))
                    
    #add to big dict
    cells_pooled_regions[brain] = pooled_regions
    
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

#get counts for all of thalamus
sois = ["Thalamus"]    

#make new dict - for all brains
all_thal_counts = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    thalamus = {}
    
    for soi in sois:
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
            thalamus[soi.name] = np.sum(np.asarray(counts))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)
            thalamus[soi] = np.sum(np.asarray(counts))
                    
    #add to big dict
    all_thal_counts[brain] = thalamus
    
#making the proper array per brain where BRAIN NAMES are removed
thal_counts = []

#initialise dummy var
i = []
for k,v in all_thal_counts.items():
    dct = all_thal_counts[k]
    for j,l in dct.items():
        i.append(l)  
    thal_counts.append(i[0])
    #re-initialise for next
    i = []    
    
#make into % counts the proper way
cell_counts_per_brain_p = np.asarray([(brain/thal_counts[i]*100) for i, brain in enumerate(cell_counts_per_brain)])
    
#%%
#pool based on tp's suggestions

#pool 1
cell_counts_per_brain_pool = np.asarray([[x[0]+x[1], x[2]+x[3]+x[4]+x[5], x[6]+x[7]+x[8], 
                                         x[9]+x[10], x[11]+x[12], 
                                         x[13]+x[14]+x[15]+x[16]+x[17]] for x in cell_counts_per_brain_p])

#%%
#pooled injections
ak_pool = np.asarray(["Lob. III, IV-V", "Lob. VIa, VIb, VII-X", 
                 "Simplex", "Crus I", "Crus II", "PM, CP"])

#pooling injection regions
expr_all_as_frac_of_lob_pool = np.asarray([[brain[0]+brain[1], brain[2]+brain[3]+brain[4]+brain[5]+brain[6]+brain[7], 
                                            brain[8], brain[9], brain[10], brain[11]+brain[12]] for brain in expr_all_as_frac_of_lob])

expr_all_as_frac_of_inj_pool = np.asarray([[brain[0]+brain[1], brain[2]+brain[3]+brain[4]+brain[5]+brain[6]+brain[7], 
                                            brain[8], brain[9], brain[10], brain[11]+brain[12]] for brain in expr_all_as_frac_of_inj])
    
primary_pool = np.asarray([np.argmax(e) for e in expr_all_as_frac_of_inj_pool])
primary_lob_n = np.asarray([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])

print(expr_all_as_frac_of_lob_pool.shape)

#normalise inj
total_lob_sum = np.asarray([np.sum(expr_all_as_frac_of_lob_pool[i]) for i in range(len(expr_all_as_frac_of_lob))])    
expr_all_as_frac_of_lob_pool_norm = np.asarray([expr_all_as_frac_of_lob_pool[i]/total_lob_sum[i] for i in range(len(expr_all_as_frac_of_lob))])

#name thalamic pools something?
#regions = np.asarray(["Parafascicular\n& Posterior complex", 
#                       "Lateral habenula,\nLateral posterior,\nPosterior triangular,\nLateral dorsal",
#                       "Paraventricular,\nCentral lateral,\nReuniens",
#                       "Ventral posteromedial\n& posterolateral",
#                       "Mediodorsal,\nVentral lateral geniculate",
#                       "Submedial, Anteroventral\nReticular, Ventral medial\nVentral anterior-lateral"])

#%%

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(cell_counts_per_brain_p[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,7))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 6
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("Mean % of thalamic counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in nuclei], fontsize="xx-small")
dst = "/home/wanglab/Desktop"
#plt.savefig(os.path.join(dst,"thal_mean_count.pdf"), bbox_inches = "tight")
#%%
#glm
X = expr_all_as_frac_of_lob_pool_norm
Y = cell_counts_per_brain_p

mat = []
pmat = []
mat_shuf = []
p_shuf = []
ars = []
rs = []
fit = []
fit_shuf = []

for itera in range(1000):
    if itera%100 == 0: print(itera)
    if itera == 0:
        shuffle = False
        inj = X.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        fit_shuf.append([])
        inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:]
    for count, region in zip(Y.T, nuclei):
        try:
            inj_ = inj[~np.isnan(count)]
            count = count[~np.isnan(count)]
    
            # intercept:
            inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
            
#            glm = sm.OLS(count, inj_)
            glm = sm.GLM(count, inj_, family=sm.families.Poisson())
            res = glm.fit()
            
            coef = res.params[:-1]
            se = res.bse[:-1]
            pvals = res.pvalues[:-1] 
    
            val = coef/se
    
            if not shuffle:
                mat.append(val)
                pmat.append(pvals)
                fit.append(res.fittedvalues)
            elif shuffle:
                mat_shuf[-1].append(val)
                p_shuf[-1].append(pvals)
                fit_shuf[-1].append(res.fittedvalues)
        except:
            inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:] #sometimes the shuffle messes stuff up
            inj_ = inj[~np.isnan(count)]
            count = count[~np.isnan(count)]
    
            # intercept:
            inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
            
#            glm = sm.OLS(count, inj_)
            glm = sm.GLM(count, inj_, family=sm.families.Poisson())
            res = glm.fit()
            
            coef = res.params[:-1]
            se = res.bse[:-1]
            pvals = res.pvalues[:-1] 
    
            val = coef/se
    
            if not shuffle:
                mat.append(val)
                pmat.append(pvals)
#                ars.append(res.rsquared_adj)
#                rs.append(res.rsquared)
            elif shuffle:
                mat_shuf[-1].append(val)
                p_shuf[-1].append(pvals)
        # inspect residuals
        """
        if not shuffle:
            plt.clf()
            plt.scatter(res.fittedvalues, res.resid)
            plt.hlines(0, res.fittedvalues.min(), res.fittedvalues.max())
            plt.title(region)
            plt.xlabel("Fit-vals")
            plt.ylabel("Residuals")
            plt.savefig(os.path.join(dst, "resid_inspection-{}.png").format(region))
        """
    

mat = np.array(mat) # region x inj
mat_shuf = np.array(mat_shuf) # region x inj
pmat = np.array(pmat) # region x inj
p_shuf = np.array(p_shuf)
fit = np.array(fit)
fit_shuf = np.array(fit_shuf)
#%%
## display
fig = plt.figure(figsize=(5,7))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat#np.median(mat_shuf, axis = 0) # NOTE abs

vmin = 0
vmax = 3
cmap = plt.cm.Reds
cmap.set_under('w')
cmap.set_over('maroon')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.2, aspect=10)
cb.set_label("Weight / SE", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# signif
sig = pmat < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05", ha="center", va="center", fontsize="small", transform=ax.transAxes)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
# xticksjt -t monokai -m 200
ax.set_xticks(np.arange(len(ak_pool))+.5)

#remaking labeles so it doesn't look squished
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in nuclei], fontsize="x-small")
dst = "/home/wanglab/Desktop"
#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")

#%%

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_fit = np.asarray([np.mean(fit.T[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,7))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_fit.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 6
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("Mean % of thalamic counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")
  
# aesthetics
# xticksjt -t monokai -m 200
ax.set_xticks(np.arange(len(ak_pool))+.5)

#remaking labeles so it doesn't look squished
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in nuclei], fontsize="x-small")
dst = "/home/wanglab/Desktop"
#plt.savefig(os.path.join(dst,"thal_mean_fit.svg"), bbox_inches = "tight")

#%%
#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_shuf_fit = np.asarray([np.mean(fit_shuf.mean(axis = 0).T[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,7))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_shuf_fit.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 6
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("Mean % of thalamic counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in nuclei], fontsize="xx-small")
dst = "/home/wanglab/Desktop"
#plt.savefig(os.path.join(dst,"thal_mean_fit_shuf.svg"), bbox_inches = "tight")

#%%

#make boxplots to make fit and shuffle fit vs. actual?
short_nuclei = ['Parafascicular n.',
 'P. complex',
 'P. triangular n.',
 'LP',
 'L. habenula',
 'LD',
 'CL',
 'PV',
 'Reuniens',
 'MD',
 'V. of lat. gen. complex',
 'VPL',
 'VPM',
 'Submedial n.',
 'RTN',
 'VM',
 'AV',
 'VA-L']

df = pd.DataFrame()
df["count"] = cell_counts_per_brain_p.T.ravel()
df["fit"] = fit.ravel()
df["shuffle"] = fit_shuf.mean(axis = 0).ravel()
df["shuffle_log"] = np.log10(fit_shuf.mean(axis = 0).ravel())
df["fit_log"] = np.log10(df["count"].values)
df["brain"] = np.array(brains*18)
df["region"] = np.repeat(np.asarray(short_nuclei), 23)
df["inj"] = np.array([ak_pool[idx] for idx in primary_pool]*18)

sns.set_style("white")

colors = ["r", "darkgoldenrod", "darkblue", "olivedrab", "darkslategray", "saddlebrown", "purple"]

g = sns.FacetGrid(df, col="inj", height=5, aspect=.5, palette = colors)

g.map(sns.barplot, "fit", "region", facecolor=(1, 1, 1, 0), ci = "sd", capsize = 0.2, edgecolor = "darkblue") 
g.map(sns.swarmplot, "fit", "region", color = "darkblue", size = 3) 
g.map(sns.barplot, "shuffle", "region", alpha = 0.4, color = "gray", ci = None) 

g.set_xlabels("% of thalamic counts")

sns.despine(offset=10)
    
plt.savefig(os.path.join(dst,"boxplot.pdf"), bbox_inches = "tight")