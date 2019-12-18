#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 19:22:25 2019

@author: wanglab
"""


from tools.analysis.network_analysis import make_structure_objects
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, pandas as pd
from skimage.external import tifffile

#custom
inj_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
cells_regions_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_contra_counts_23_brains.csv"
dst = "/home/wanglab/Desktop"

#making dictionary of injection sites
injections = {}

#MAKE SURE THEY ARE IN THIS ORDER
brains = ["20170410_tp_bl6_lob6a_ml_repro_01",
#         "20160823_tp_bl6_cri_500r_02",
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
expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()] 
                        for mouse in inj_raw])    
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for nm, lob in atlas_rois.items()] 
                        for mouse in inj_raw])    
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

print(expr_all_as_frac_of_lob[15])

#%%
#making dictionary of cells by region
cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
#pooling regions
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#%%

#GET ONLY THALAMIC POOLS
nuclei = ["Thalamus",
     "Ventral posteromedial nucleus of the thalamus",
     "Ventral posterolateral nucleus of the thalamus",
     "Ventral anterior-lateral complex of the thalamus",
     "Anteroventral nucleus of thalamus",
     "Lateral dorsal nucleus of thalamus",
     "Paraventricular nucleus of the thalamus",
     "Medial habenula",
     "Lateral posterior nucleus of the thalamus",
     "Posterior triangular thalamic nucleus",
     "Mediodorsal nucleus of thalamus",
     "Posterior complex of the thalamus",
     "Ventral medial nucleus of the thalamus",
     "Reticular nucleus of the thalamus"
]

#make new dict - for all brains
cells_pooled_regions = {} #for raw counts

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    
    for soi in nuclei:
        counts = []
        try:
            soi = [s for s in structures if s.name==soi][0]
            counts.append(cells_regions.loc[cells_regions.Structure == soi.name, brain].values[0]) #store counts in this list
            #add to volume list from LUT
            progeny = [str(xx.name) for xx in soi.progeny]
            #now sum up progeny
            if len(progeny) > 0:
                for progen in progeny:
                    counts.append(cells_regions.loc[cells_regions.Structure == progen, brain].values[0])
                            #add to volume list from LUT
            c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
        except:
            counts.append(cells_regions.loc[cells_regions.Structure == soi, brain].values[0])                    
            c_pooled_regions[soi] = np.sum(np.asarray(counts))
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions

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
#make into % counts the proper way
pcounts = np.asarray([((brain[1:]/brain[0])*100) for brain in cell_counts_per_brain])    
    
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
expr_all_as_frac_of_inj_pool_norm = np.asarray([brain/brain.sum() for brain in expr_all_as_frac_of_inj_pool])

#change vars (consistent with NC)
regions = np.asarray(nuclei)[1:]
#%%
#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(pcounts[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 8
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,(vmax-vmin)/2 +1)
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
## exact value annotations
#for ri,row in enumerate(show):
#    for ci,col in enumerate(row):
#        pass
#        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(nuclei[1:]))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in nuclei[1:]], fontsize="small")
dst = "/home/wanglab/Desktop"
#plt.savefig(os.path.join(dst,"thal_mean_count.pdf"), bbox_inches = "tight")

#%%
#glm
X = expr_all_as_frac_of_inj_pool_norm
Y = pcounts

c_mat = []
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
                c_mat.append(coef)
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
fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat

vmin = 0
vmax = 5
whitetext = 4
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
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
        if col > whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")

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
ax.set_xticks(np.arange(len(ak_pool))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")

plt.savefig(os.path.join(dst,"thal_glm.pdf"), bbox_inches = "tight")