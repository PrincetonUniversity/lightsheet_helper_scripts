#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:10:56 2019

@author: wanglab
"""

from scipy.stats import spearmanr
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl
from skimage.external import tifffile
from scipy import ndimage
import unicodedata
import string

valid_filename_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
char_limit = 255

def clean_filename(filename, whitelist=valid_filename_chars, replace=' '):
    """ from https://gist.github.com/wassname/1393c4a57cfcbf03641dbc31886123b8 """
    # replace spaces
    for r in replace:
        filename = filename.replace(r,'_')
    
    # keep only valid ascii chars
    cleaned_filename = unicodedata.normalize('NFKD', filename).encode('ASCII', 'ignore').decode()
    
    # keep only whitelisted chars
    cleaned_filename = ''.join(c for c in cleaned_filename if c in whitelist)
    if len(cleaned_filename)>char_limit:
        print("Warning, filename truncated because it was over {}. Filenames may no longer be unique".format(char_limit))
        return cleaned_filename[:char_limit] 

#custom
inj_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/injection_sites"
atl_pth = "/jukebox/LightSheetTransfer/atlas/cb_sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/cb_annotation_sagittal_atlas_20um_iso.tif"
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
    
#make it an array   
inj_raw = np.array([inj.astype(bool) for nm, inj in injections.items()])

#making dictionary of all 3d cells
cells_coords = {}

for pth in os.listdir(cells_pth):
    print(pth)
    pth = os.path.join(cells_pth, pth)
    cells = np.load(pth)
    print(cells.shape)
    cells_coords[os.path.basename(pth)[:-4]] = [str(tuple(cell)) for cell in cells]

#make it an array
cell_coords_raw = np.array([coords for nm, coords in cells_coords.items()])

#from badura et al.
##  glm
mat = []
pmat = []
mat_shuf = []
p_shuf = []
ars = []
for itera in range(100):
    print(itera)
    if itera == 0:
        shuffle = False
        inj = inj_raw.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        inj = inj_raw[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:]
    for i, coord in enumerate(cell_coords_raw):
        inj_ = inj_raw[i]
        coord = coord[i]

        # intercept:
        inj_ = np.concatenate([inj_, np.ones(inj_.shape)], axis=1)

        glm = sm.OLS(coord, inj_)
        #glm_fam = sm.families.Gaussian(sm.genmod.families.links.identity)
        #glm = sm.GLM(b, ex_, family=glm_fam)
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
        """
        if not shuffle:
            pl.clf()
            pl.scatter(res.fittedvalues, res.resid)
            pl.hlines(0, res.fittedvalues.min(), res.fittedvalues.max())
            pl.title(bname)
            pl.xlabel("Fit-vals")
            pl.ylabel("Residuals")
            pl.savefig("figs/resid_inspection-{}.png".format(bname))
        """

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
reg = regions[ordr]

## display
fig = plt.figure(figsize=(4,7))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = np.abs(mat) # NOTE abs
amax = 4.47 #computed from juvenile, so colors match for both #np.max(np.abs(mat))
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
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", fontsize=8, ha="center", va="center")

# signif
sig = pmat < .05#/np.size(pmat)
null = (p_shuf < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    #ax.text(x+.5, y+.5, "x", fontsize=15, ha="center", va="center", transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) X"s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="small", transform=ax.transAxes)

# aesthetics
# xticks
ax.set_xticks(np.arange(len(ak))+.5)
ax.set_xticklabels(ak, rotation=30, fontsize="x-small", ha="center")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,bn)], fontsize="x-small")
ax.set_yticklabels(["{}".format(bi) for ar,bi in zip(ars, cell_counts_per_brain)], fontsize="x-small")