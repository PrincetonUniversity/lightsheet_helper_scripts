#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 20:09:09 2019

@author: wanglab
"""
from lightsheet.network_analysis import make_structure_objects
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl
import tifffile

#VARIABLES FOR GLM   
#cell_counts_per_brain = np.load("/jukebox/wang/zahra/modeling/h129/thalamus/cell_counts_per_brain.npy")
#
##remove non associated areas
#mask = np.ones(cell_counts_per_brain[0].shape, dtype=bool)
#mask[[0, 1, 5, 15, 16]] = False
#cell_counts_per_brain = np.asarray([xx[mask] for xx in cell_counts_per_brain])


#group regions differently
cells_regions_pth = '/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/cell_count_by_region/dataframe.p'
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_200um_edges_only.tif"

cells_regions = pckl.load(open(cells_regions_pth, "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

brains = list(cells_regions.keys())

###########################################################WORKING ON THIS###########################################################
#pooling regions
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#GET ONLY NEOCORTICAL POOLS
sois = ["Thalamus, sensory-motor related", "Thalamus, polymodal association cortex related"]

#make new dict - for all brains
cells_pooled_regions = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    pooled_regions = {}
    
    for soi in sois:
        try:
            soi = [s for s in structures if s.name==soi][0]
            print(soi.name)
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
    
injp = np.load("/jukebox/wang/zahra/modeling/h129/thalamus/fraction_of_lobule_pc_layer.npy")

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
#LOOKING AT PERCENT CELL COUNTS
cell_counts_per_brain = np.asarray([(xx/np.sum(xx))*100 for xx in cell_counts_per_brain])
#np.save("/jukebox/wang/zahra/modeling/h129/thalamus/pooled_cell_counts_per_brain.npy", cell_counts_per_brain)

#pooling injection regions
#expr_all_as_frac_of_lob_pool = []
#for brain in injp:
#    #group 1: ant vermis; group 2: post vermis; group 3: post hemisphere
#    pool_brain = [brain[0], brain[1]+brain[2]+brain[3]+brain[4]+brain[5]+brain[6], 
#                  brain[7]+brain[8]+brain[9]+brain[10]+brain[11]]
#    
#    expr_all_as_frac_of_lob_pool.append(pool_brain)    
#
#expr_all_as_frac_of_lob_pool = np.asarray(expr_all_as_frac_of_lob_pool)

#regions = np.load("/jukebox/wang/zahra/modeling/h129/thalamus/thalamic_regions.npy")
regions = sois
##remove non associated areas
#mask = np.ones(regions.shape[0], dtype=bool)
#mask[[0, 1, 5, 15, 16]] = False
#regions = regions[mask]

ak = np.load("/jukebox/wang/zahra/modeling/h129/thalamus/cerebellar_regions.npy")
#ak = np.asarray(["Anterior vermis", "Posterior vermis", "Posterior hemisphere"])
#injp = expr_all_as_frac_of_lob_pool

#np.random.shuffle(cell_counts_per_brain)
#%%
#from badura et al.
##  glm
mat = []
pmat = []
mat_shuf = []
p_shuf = []
ars = []
for itera in range(1000):
    print(itera)
    if itera == 0:
        shuffle = False
        inj = injp.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        inj = injp[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:]
    for count, region in zip(cell_counts_per_brain.T, regions):
        inj_ = inj[~np.isnan(count)]
        count = count[~np.isnan(count)]

        # intercept:
        inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]], axis=1)

        glm = sm.OLS(count, inj_)
#        glm_fam = sm.families.Gaussian(sm.genmod.families.links.identity)
#        glm = sm.GLM(count, inj_, family=glm_fam)
        res = glm.fit()

        coef = res.params[:-1] 
        se = res.bse[:-1]
        pvals = res.pvalues[:-1]

        val = coef/se

        if not shuffle:
            mat.append(val)
            pmat.append(pvals)
            ars.append(res.rsquared_adj)
        elif shuffle:
            mat_shuf[-1].append(val)
            p_shuf[-1].append(pvals)
        
        # inspect residuals
#        if not shuffle:
#            plt.clf()
#            plt.scatter(res.fittedvalues, res.resid)
#            plt.hlines(0, res.fittedvalues.min(), res.fittedvalues.max())
#            plt.title(region)
#            plt.xlabel("Fit-vals")
#            plt.ylabel("Residuals")
#            plt.savefig("/home/wanglab/Desktop/resid_inspection-{}.png".format(region))
        

mat = np.array(mat) # region x inj
mat_shuf = np.array(mat_shuf) # region x inj
pmat = np.array(pmat) # region x inj
p_shuf = np.array(p_shuf)
ars = np.array(ars)

#ordr = np.argsort(np.max(np.abs(mat),axis=1))
ordr = np.arange(len(mat))
mat = mat[ordr,:]
mat_shuf = mat_shuf[:,ordr,:]
p_shuf = p_shuf[:,ordr,:]
pmat = pmat[ordr,:]
reg = regions[ordr]

#%%
## display
fig = plt.figure()#figsize=(4,7))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat # NOTE abs
amax = 4.7 #computed from np.max(np.abs(mat))
#amax = np.max(np.abs(mat))
#vmin = -amax
vmin = np.min(mat)-0.5
vmax = np.max(mat)+0.5
cmap = plt.cm.RdBu_r
#cmap = plt.cm.coolwarm

# discrete colorbar details
bounds = np.linspace(-5,5,11)
#bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", shrink=0.5, aspect=10)
cb.set_label("| Weight / SE |", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# signif
sig = pmat < .05#/np.size(pmat)
null = (p_shuf < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x+.5, y+.5, "x", fontsize=15, ha="center", va="center", transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05".format(nullmean, nullstd), ha="center", va="center", fontsize="small", transform=ax.transAxes)
ax.text(.5, 1.06, "X: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) X's are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
# xticks
ax.set_xticks(np.arange(len(ak))+.5)

#remaking labeles so it doesn't look squished
lbls = ['Lob. IV-V',
       'Lob. VIa', 'Lob. VIb', 'Lob. VII', 'Lob. VIII',
       'Lob, IX', 'Lobule X', 'Simplex', 'Crus 1', 'Crus 2',
       'Paramedian', 'Copula pyramidis']
lbls = np.asarray(lbls)
ax.set_xticklabels(lbls, rotation=30, fontsize=4, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,bn)], fontsize="x-small")
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="x-small")
plt.savefig("/home/wanglab/Desktop/weights.pdf")