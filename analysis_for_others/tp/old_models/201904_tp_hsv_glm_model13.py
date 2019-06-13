#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 17:55:33 2019

@author: wanglab
"""

from lightsheet.network_analysis import make_structure_objects
from scipy.stats import spearmanr
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl
from skimage.external import tifffile

#custom
inj_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
#making dictionary of injection sites
injections = {}

for pth in os.listdir(inj_pth):
    print(pth)
    pth = os.path.join(inj_pth, pth)
    injection = tifffile.imread(pth)
    print(injection.shape)
    injections[os.path.basename(pth)[:-8]] = injection #files have 2 .tif in the end
    
inj_raw = np.array([inj.astype(bool) for nm, inj in injections.items()])
    
atl_raw = tifffile.imread(atl_pth)
ann_raw = tifffile.imread(ann_pth)[:, 423:, :] #apparent cutoff
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)

#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lobule IV-V": 1091,
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

frac2 = np.array([ms[scd] for ms,scd in zip(expr_all_as_frac_of_lob, secondary)])
print("\nFraction of lobule for secondary lobules: {:0.2f} ± {:0.2f}".format(frac2.mean(), frac2.std(ddof=1)))

print("\nGropuing by primary lobule, mean fraction of secondary lobule covered:")
for p in np.unique(primary):
    scd = secondary[primary==p]
    mc = expr_all_as_frac_of_lob[primary==p]
    fracs = [mi[si] for mi,si in zip(mc,scd)]
    print("primary =",list(iids.keys())[p])
    print("\tfraction of secondary = {:0.2f}".format(np.mean(fracs)))


#save out
np.save("/jukebox/wang/zahra/modeling/h129/thalamus/fraction_of_lobule.npy", expr_all_as_frac_of_lob)
#%%
#VARIABLES FOR GLM   
cell_counts_per_brain = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/cell_counts_per_brain.npy")
#LOOKING AT PERCENT CELL COUNTS
cell_counts_per_brain = np.asarray([(xx/np.sum(xx))*100 for xx in cell_counts_per_brain])

injp = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/fraction_of_lobule_pc_layer.npy")

#pooling injection regions
expr_all_as_frac_of_lob_pool = []
for brain in injp:
    #group 1: ant vermis; group 2: post vermis; group 3: post hemisphere
    pool_brain = [brain[0], brain[1]+brain[2]+brain[3]+brain[4]+brain[5]+brain[6], 
                  brain[7]+brain[8]+brain[9]+brain[10]+brain[11]]
    
    expr_all_as_frac_of_lob_pool.append(pool_brain)    

expr_all_as_frac_of_lob_pool = np.asarray(expr_all_as_frac_of_lob_pool)  

regions = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/neocortical_regions.npy")

##remove non associated areas
#mask = np.ones(regions.shape[0], dtype=bool)
#mask[[0, 1, 5, 15, 16]] = False
#regions = regions[mask]

#ak = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/cerebellar_regions.npy")
ak = np.asarray(["Ant. vermis", "Post. vermis", "Post. hemisphere"])

#np.random.shuffle(cell_counts_per_brain)
injp = expr_all_as_frac_of_lob_pool

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
#amax = 0.59 #computed from np.max(np.abs(mat))
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
#
##remaking labeles so it doesn't look squished
#lbls = ['Lob. IV-V',
#       'Lob. VIa', 'Lob. VIb', 'Lob. VII', 'Lob. VIII',
#       'Lob, IX', 'Lobule X', 'Simplex', 'Crus 1', 'Crus 2',
#       'Paramedian', 'Copula pyramidis']
#lbls = np.asarray(lbls)
ax.set_xticklabels(ak, rotation=20, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,bn)], fontsize="x-small")
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="x-small")
plt.savefig("/home/wanglab/Desktop/weights.svg")

#%%
#plot r squared adj
np.save("/jukebox/wang/zahra/modeling/h129/neocortex/r2adjusted_pc_layer_pooled_inj.npy", ars)

fig, ax = plt.subplots(figsize=(11,8))

ax.scatter(range(17), ars, color = "black")
ax.axhline(0, color = 'grey')
ax.set_xticks(np.arange(len(regions)))
ax.set_xticklabels(regions, rotation=20, fontsize=8, ha="right")
ax.set_xlabel("Thalamic regions")
ax.set_ylabel("Adjusted R-squared")
plt.savefig('/home/wanglab/Desktop/r2sdj.svg')