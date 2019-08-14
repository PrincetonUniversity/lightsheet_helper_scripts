#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 12:28:11 2019

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
inj_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/injection_sites"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_200um_edges_only.tif"
cells_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/cell_count_by_coordinate_only_including_structure"
cells_regions_pth = '/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/cell_count_by_region/nc_dataframe.p'

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

## compute stat
for lid,lstr in zip([3, 4, 5, 7, 8, 10, 11, 12, 13, 14],["Lobule IV-V", "Lobule VIa", "Lobule VIb","Lobule VIII", "Lobule IX", "Simplex lobule", "Crus 1", "Crus 2", "Paramedian lobule", "Copula pyramidis"]):
    grp = expr_all_as_frac_of_lob[primary==lid][:,lid]
    mean = grp.mean()
    std = grp.std(ddof=1)
    print("{}\t{:0.2f} ± {:0.2f}".format(lstr, mean, std))

frac2 = np.array([ms[scd] for ms,scd in zip(expr_all_as_frac_of_lob, secondary)])
print("\nFraction of lobule for secondary lobules: {:0.2f} ± {:0.2f}".format(frac2.mean(), frac2.std(ddof=1)))

print("\nGropuing by primary lobule, mean fraction of secondary lobule covered:")
for p in np.unique(primary):
    scd = secondary[primary==p]
    mc = expr_all_as_frac_of_lob[primary==p]
    fracs = [mi[si] for mi,si in zip(mc,scd)]
    print("primary =",list(iids.keys())[p])
    print("\tfraction of secondary = {:0.2f}".format(np.mean(fracs)))

#making dictionary of all 3d cells
cells_coords = {}

for pth in os.listdir(cells_pth):
    print(pth)
    pth = os.path.join(cells_pth, pth)
    cell = np.load(pth)
    print(cell.shape)
    cells_coords[os.path.basename(pth)[:-4]] = cell #files have 2 .tif in the end

#making dictionary of cells by region
cells_regions = pckl.load(open(cells_regions_pth, "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

brains = list(cells_regions.keys())

###########################################################WORKING ON THIS###########################################################
#pooling regions
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#GET ONLY NEOCORTICAL POOLS
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#make new dict - for all brains
cells_pooled_regions = {}

for brain in cells_regions.keys():    
    #make new dict - this is for EACH BRAIN
    pooled_regions = {}
    
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print(soi.name)
        progeny = [str(xx.name) for xx in soi.progeny]
        print(progeny)
        counts = [] #store counts in this list
        if len(progeny) > 0:
            for progen in progeny:
                for k, v in cells_regions[brain].items():
                    if k == progen:
                        counts.append(v)
        pooled_regions[soi.name] = np.sum(np.asarray(counts))
                    
    #add to big dict
    cells_pooled_regions[brain] = pooled_regions
    
###########################################################WORKING ON THIS###########################################################

#making the proper array per brain where regions are removed
cell_counts_per_brain = []

#initialise dummy var
i = []
for k,v in cells_pooled_regions.items():
    dct = cells_pooled_regions[k]
    for j,l in dct.items():
      i.append(l)  
    cell_counts_per_brain.append(i)
    #re-initialise for nect
    i = []

#%%
#VARIABLES FOR GLM   
cell_counts_per_brain = np.asarray(cell_counts_per_brain)        
        
injp = np.nan_to_num(expr_all_as_frac_of_lob)

regions = np.asarray(sois)

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
        inj_ = injp[~np.isnan(count)]
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
        

mat = np.array(mat) # behavior x voxel
mat_shuf = np.array(mat_shuf) # behavior x voxel
pmat = np.array(pmat) # behavior x voxel
p_shuf = np.array(p_shuf)
ars = np.array(ars)

#ordr = np.argsort(np.max(np.abs(mat),axis=1))
ordr = np.arange(len(mat))
mat = mat[ordr,:]
mat_shuf = mat_shuf[:,ordr,:]
p_shuf = p_shuf[:,ordr,:]
pmat = pmat[ordr,:]
reg = regions

#%%
## display
fig = plt.figure()#figsize=(4,7))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = np.abs(mat) # NOTE abs
amax = 4.7 #computed from np.max(np.abs(mat))
#amax = np.max(np.abs(mat))
#vmin = -amax
vmin = 0
vmax = amax
#cmap = pl.cm.RdBu_r
cmap = plt.cm.Reds

# discrete colorbar details
#bounds = np.linspace(-5,5,11)
bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", shrink=0.5, aspect=10)
cb.set_label("| Weight / SE |", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)

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
lbls = ['Lob. I', 'Lob. II', 'Lob. III', 'Lob. IV-V',
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
#%%
## specific relevant scatter plots

rchoices = regions.copy()

echoices = np.asarray(['Lobule IV-V','Lobule VIa','Lobule VIII','Lobule IX', 'Crus 1', 'Crus 2', 'Paramedian lobule'])
primids = [3, 4, 7, 8, 11, 12, 13]
prim = primary.copy()

choice_idx = 1 # from list directly above this, default 0
echoice = echoices[choice_idx]
col = (prim==primids[choice_idx]).astype(int) # should match echoice

fig,ax = plt.subplots(1,1,figsize=(5,5), gridspec_kw=dict(bottom=.2, left=.2, right=.9, top=.9))
for rg in rchoices:
    ax.clear()

    bi = int(np.where(regions == rg)[0].squeeze())
    ei = int(np.where(np.array(ak) == echoice)[0].squeeze())
    ex_vec = inj_[:,ei]
    b_vec = cell_counts_per_brain[:,bi]
    if not np.nonzero(b_vec)[0].shape == (0,):
        mfc = ["none","k"]
    
        for mfci,idx in zip(mfc, [0,1]):
            ax.scatter(ex_vec[col==idx], b_vec[col==idx], color="k", facecolor=mfci, linewidth=.2, s=160)
            ax.set_xlabel('Fraction of {} Molecular layer labelled'.format(echoice))
            ax.set_ylabel("{}".format(rg))
        plt.savefig('/home/wanglab/Desktop/figs/scatter/scatter_{}_{}.pdf'.format(rg, echoice))

## save out all inspection scatterplots
plt.close("all")

for rname, reh in zip(regions, cell_counts_per_brain.T):
    if not np.nonzero(reh)[0].shape == (0,):
        print(rname)
        fig,axs = plt.subplots(3,4)
        for lob,lobname,ax in zip(inj.T, ak, axs.ravel()):
            ax.scatter(lob, reh)
            ax.set_title(lobname, fontsize="xx-small")
            ax.set_xticks([])
            ax.set_yticks([])
        plt.figtext(.5, .05, "Fraction of Molecular layer labelled", fontsize="small")
        plt.figtext(.1, .5, rname, fontsize="small", rotation=90)
        fig.savefig("/home/wanglab/Desktop/figs/byreg/byreg-"+rname+".pdf")
        plt.close(fig)

#%%
## correlation matrices
lname1 = "Lobule VIa"
lname2 = "Crus 1"
rho0 = spearmanr(cell_counts_per_brain[(ak[primary]==lname1)], cell_counts_per_brain[(ak[primary]==lname2)])[0]
rho0[np.isnan(rho0)] = 0
fig,ax = plt.subplots(1,1)
ax.pcolor(np.abs(rho0))
ax.set_title("Correlation matrices", fontsize="small")