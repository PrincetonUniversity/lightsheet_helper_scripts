#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 19:06:33 2019

@author: wanglab

outputs all figures and dataframes to 'dst' directory
"""

import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, matplotlib as mpl, pandas as pd
from scipy.stats import median_absolute_deviation as mad

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#import data
main_data_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_contra_ipsi_counts_densities.p"
data = pckl.load(open(main_data_pth, "rb"), encoding = "latin1")

#injection site analysis
inj_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data.p"
inj_dct = pckl.load(open(inj_pth, "rb"), encoding = "latin1")

#inj volumes
inj_vol_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_inj_vol.p" 
inj_vol_dct = pckl.load(open(inj_vol_pth, "rb"), encoding = "latin1")
inj_vol = inj_vol_dct["thal_inj_vol"]
iv = []
for k,v in inj_vol.items():
    iv.append(v)
vols = [xx/1e4 for xx in iv]

#set dst 
sv_dst = "/home/wanglab/Desktop"

#mask unwanted brains - dropping fast spreading ones and a low count one
curated_brains = [True, False, False, True, True, True, False, True, True, True, True, True,
                   True, True, True, True, True, True, True, True, True, True, True]

cell_counts_per_brain_left = data["cell_counts_per_brain_left"][curated_brains]
cell_counts_per_brain_right = data["cell_counts_per_brain_right"][curated_brains]
density_per_brain_left = data["density_per_brain_left"][curated_brains]
density_per_brain_right = data["density_per_brain_right"][curated_brains]
volume_per_brain = data["volume_per_brain_left"][curated_brains]
lr_dist = data["lr_dist"]

brains = np.array(inj_dct["brainnames"])[curated_brains]
primary_pool = inj_dct["primary_pool"][curated_brains]
ak_pool = inj_dct["cb_regions_pool"]
inj = inj_dct["expr_all_as_frac_of_inj_pool"][curated_brains]
vols = np.array(vols)[curated_brains]

#-------------------------------------------------------------------------------------------------------------------------------------
#preprocessing
scale_factor = 0.020
thal_left_counts = cell_counts_per_brain_left
thal_right_counts = cell_counts_per_brain_right
thal_density_left = density_per_brain_left
thal_density_right = density_per_brain_right

lrv = np.array(list(lr_dist.values()))[curated_brains]
lr_brains = np.array(list(lr_dist.keys()))[curated_brains]

_ccontra = []; _cipsi = []; _dcontra = []; _dipsi = []
for i in range(len(lr_brains)):
    if lrv[i] > 0: #right
        #counts
        _ccontra.append(thal_left_counts[i])
        _cipsi.append(thal_right_counts[i])
        #density
        _dcontra.append(thal_density_left[i])
        _dipsi.append(thal_density_right[i])
    elif lrv[i] < 0: #left
        #counts
        _ccontra.append(thal_right_counts[i])
        _cipsi.append(thal_left_counts[i])
        #density
        _dcontra.append(thal_density_right[i])
        _dipsi.append(thal_density_left[i])


_ccontra = np.asarray(_ccontra).T; _dcontra = np.asarray(_dcontra).T
_cipsi = np.asarray(_cipsi).T; _dipsi = np.asarray(_dipsi).T
_dratio = np.asarray([_dcontra[i]/_dipsi[i] for i in range(len(_dcontra))])
_cratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
#make into one
_dist = lrv
 
_inj = np.asarray([inj[i] for i in range(len(inj)) if brains[i] in lr_brains])
_primary_pool = np.asarray([primary_pool[i] for i in range(len(primary_pool)) if brains[i] in lr_brains])

#sort by injection fractions
sort_primary = np.sort(primary_pool)
sort_var = np.argsort(primary_pool, axis = 0)
sort_inj = _inj[sort_var]
sort_dist = _dist[sort_var]
sort_ccontra = _ccontra.T[sort_var]
sort_cipsi = _cipsi.T[sort_var]
sort_cratio = _cratio.T[sort_var]
sort_dcontra = _dcontra.T[sort_var]
sort_dipsi = _dipsi.T[sort_var]
sort_dratio = _dratio.T[sort_var]
sort_vox_per_region = volume_per_brain[np.argsort(_dist, axis = 0)]

sort_brains = np.array(lr_brains)[np.argsort(primary_pool)]
sort_vols = np.array(vols)[np.argsort(primary_pool)]

print(sort_dist.shape)
print(sort_cratio.shape)

#group thalamus regions in smaller, meta regions
grps = np.array(["Sensory-motor thalamus" , "Polymodal thalamus"])
sort_ccontra_pool = np.asarray([[np.sum(xx[:5]), np.sum(xx[5:])] for xx in sort_ccontra])
sort_dcontra_pool = np.asarray([[np.sum(xx[:5]), np.sum(xx[5:])] for xx in sort_ccontra])/(np.asarray([[np.sum(xx[:5]), 
                                 np.sum(xx[5:])] for xx in sort_vox_per_region])*(scale_factor**3))
sort_cipsi_pool = np.asarray([[np.sum(xx[:5]), np.sum(xx[5:])] for xx in sort_cipsi])
sort_dipsi_pool = np.asarray([[np.sum(xx[:5]), np.sum(xx[5:])] for xx in sort_cipsi])/(np.asarray([[np.sum(xx[:5]), 
                                 np.sum(xx[5:])] for xx in sort_vox_per_region])*(scale_factor**3))
sort_cratio_pool = np.asarray([sort_ccontra_pool[i]/sort_cipsi_pool[i] for i in range(len(sort_ccontra_pool))])
sort_dratio_pool = np.asarray([sort_dcontra_pool[i]/sort_dipsi_pool[i] for i in range(len(sort_dcontra_pool))])

#-------------------------------------------------------------------------------------------------------------------------------------

df = pd.DataFrame()
df["mean_dratio_per_inj_smthal"] = np.asarray([np.mean(sort_dratio_pool.T[0][np.where(_primary_pool == idx)[0]], axis=0) for 
                                         idx in np.unique(_primary_pool)])
df["mean_dratio_per_inj_polythal"] = np.asarray([np.mean(sort_dratio_pool.T[1][np.where(_primary_pool == idx)[0]], axis=0) for 
                                         idx in np.unique(_primary_pool)])
df["med_dratio_per_inj_smthal"] = np.asarray([np.median(sort_dratio_pool.T[0][np.where(_primary_pool == idx)[0]], axis=0) for 
                                         idx in np.unique(_primary_pool)])
df["med_dratio_per_inj_polythal"] = np.asarray([np.median(sort_dratio_pool.T[1][np.where(_primary_pool == idx)[0]], axis=0) for 
                                         idx in np.unique(_primary_pool)])
df["std_dratio_per_inj_smthal"] = np.asarray([np.std(sort_dratio_pool.T[0][np.where(_primary_pool == idx)[0]], axis=0) for 
                                         idx in np.unique(_primary_pool)])
df["std_dratio_per_inj_polythal"] = np.asarray([np.std(sort_dratio_pool.T[1][np.where(_primary_pool == idx)[0]], axis=0) for 
                                         idx in np.unique(_primary_pool)])
df["est_std_dratio_per_inj_smthal"] = np.asarray([mad(sort_dratio_pool.T[0][np.where(_primary_pool == idx)[0]], axis=0)/0.6745 for 
                                         idx in np.unique(_primary_pool)])
df["est_std_dratio_per_inj_polythal"] = np.asarray([mad(sort_dratio_pool.T[1][np.where(_primary_pool == idx)[0]], axis=0)/0.6745 for 
                                         idx in np.unique(_primary_pool)])    
df.index = ak_pool
df = df.round(2)
df.to_csv(os.path.join(sv_dst, "thal_contra_ipsi_ratios_by_inj.csv"))

#-------------------------------------------------------------------------------------------------------------------------------------

df = pd.DataFrame()
df["mean_dratio_per_inj_vermis"] = [np.mean(sort_dratio_pool.T[i][np.where((_primary_pool == 0) | (_primary_pool == 1))], 
                                                axis=0) for i in range(len(sort_dratio_pool.T))] #only vermis    
df["mean_dratio_per_inj_hem"] = [np.mean(sort_dratio_pool.T[i][np.where((_primary_pool != 0) & (_primary_pool != 1))], 
                                                axis=0) for i in range(len(sort_dratio_pool.T))] #only hem
df["med_dratio_per_inj_vermis"] = [np.median(sort_dratio_pool.T[i][np.where((_primary_pool == 0) | (_primary_pool == 1))], 
                                                axis=0) for i in range(len(sort_dratio_pool.T))] #only vermis    
df["med_dratio_per_inj_hem"] = [np.median(sort_dratio_pool.T[i][np.where((_primary_pool != 0) & (_primary_pool != 1))], 
                                                axis=0) for i in range(len(sort_dratio_pool.T))] #only hem
df["std_dratio_per_inj_vermis"] = [np.std(sort_dratio_pool.T[i][np.where((_primary_pool == 0) | (_primary_pool == 1))], 
                                                axis=0) for i in range(len(sort_dratio_pool.T))] #only vermis    
df["std_dratio_per_inj_hem"] = [np.std(sort_dratio_pool.T[i][np.where((_primary_pool != 0) & (_primary_pool != 1))], 
                                                axis=0) for i in range(len(sort_dratio_pool.T))] #only hem
    
df["est_std_dratio_per_inj_vermis"] = [mad(sort_dratio_pool.T[i][np.where((_primary_pool == 0) | (_primary_pool == 1))], 
                                                axis=0)/0.6745 for i in range(len(sort_dratio_pool.T))] #only vermis
df["est_std_dratio_per_inj_hem"] = [mad(sort_dratio_pool.T[i][np.where((_primary_pool != 0) & (_primary_pool != 1))], 
                                                axis=0)/0.6745 for i in range(len(sort_dratio_pool.T))] #only hem  
df.index = ["SM", "Poly"]
df = df.round(2)
df.to_csv(os.path.join(sv_dst, "thal_contra_ipsi_ratios_by_vermis_hemisphere.csv"))

#-------------------------------------------------------------------------------------------------------------------------------------
## display
fig, axes = plt.subplots(ncols = 1, nrows = 6, figsize = (10,4), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,0.5,0.8,0.8,0.8,0.5]})


#set colormap specs
vmaxcounts = 50
whitetext = 10

#inj fractions
ax = axes[0]

show = np.fliplr(sort_inj).T

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")

#inj vols
ax = axes[1]

show = np.asarray([sort_vols])

vmin = 0
vmax = 12
cmap = plt.cm.Greens 
cmap.set_over('darkgreen')
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, 
                  boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 10:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", 
                    ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", 
                    ha="center", va="center", fontsize="x-small")

ylabel = ["Injection volume\n($10^4$ voxels)"]
ax.set_yticks(np.arange(len(ylabel))+.5)
ax.set_yticklabels(ylabel, fontsize="x-small")

ax = axes[2]
show = sort_dcontra_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcounts
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=1.5, aspect=10)
cb.set_label("Density (Cells/$mm^3$)", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(-0.25,0.5)

ax = axes[3]
show = sort_dipsi_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcounts
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=1.5, aspect=10)
cb.set_label("Density (Cells/$mm^3$)", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.25,0.5)


ax = axes[4]
show = sort_dratio_pool.T
yaxis = grps

vmin = 0.7
vmax = 1.5
cmap = plt.cm.Blues
cmap.set_over("navy")
#colormap
bounds = np.linspace(vmin,vmax,5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
cb.set_label("Density ratio", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col > 1.5:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")
ax.set_ylabel("Contra/Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.25,0.5)

ax = axes[5]
show = np.asarray([sort_dist])
br = lr_brains 

vmin = -100
vmax = 80
cmap = plt.cm.RdBu_r
cmap.set_over('maroon')
cmap.set_under('midnightblue')
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=2, aspect=10)
cb.set_label("Left to right", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < -75 or col > 70:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")        

#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(sort_brains))+.5)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=5, ha="right")

ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

plt.savefig(os.path.join(sv_dst, "thal_density_ratios.pdf"), bbox_inches = "tight")


#-------------------------------------------------------------------------------------------------------------------------------------
## display
fig, axes = plt.subplots(ncols = 1, nrows = 6, figsize = (10,4), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,0.4,0.8,0.8,0.8,0.5]})

#set colormap specs
vmaxcounts = 150
whitetext = 30
    
#inj fractions
ax = axes[0]

show = np.fliplr(sort_inj).T

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")

#inj vols
ax = axes[1]

show = np.asarray([sort_vols])

vmin = 0
vmax = 8
cmap = plt.cm.Greens 
cmap.set_over('darkgreen')
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, 
                  boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 6:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", 
                    ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", 
                    ha="center", va="center", fontsize="x-small")

ylabel = ["Injection volume\n($10^4$ voxels)"]
ax.set_yticks(np.arange(len(ylabel))+.5)
ax.set_yticklabels(ylabel, fontsize="x-small")

ax = axes[2]
show = sort_ccontra_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcounts
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=1.5, aspect=10)
cb.set_label("Cell count", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(-0.25,0.5)

ax = axes[3]
show = sort_cipsi_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcounts
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.8, aspect=10)
cb.set_label("Cell count", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.25,0.5)


ax = axes[4]
show = sort_cratio_pool.T
yaxis = grps

vmin = 0.7
vmax = 1.5
cmap = plt.cm.Blues
cmap.set_over("navy")
#colormap
bounds = np.linspace(vmin,vmax,5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
cb.set_label("Cell count ratio", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col > 1.5:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra/Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.25,0.5)

ax = axes[5]
show = np.asarray([sort_dist])
br = lr_brains 

vmin = -100
vmax = 80
cmap = plt.cm.RdBu_r
cmap.set_over('maroon')
cmap.set_under('midnightblue')
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=2, aspect=10)
cb.set_label("Left to right", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < -75 or col > 70:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")        

#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(br))+.5)
lbls = np.asarray(br)
ax.set_xticklabels(br, rotation=30, fontsize=5, ha="right")

ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

plt.savefig(os.path.join(sv_dst, "thal_count_ratios.pdf"), bbox_inches = "tight")
#-------------------------------------------------------------------------------------------------------------------------------------
#basic statistics for these ratios

df = pd.DataFrame()
#decimal to round by
d = 2
df["median_density_ratio"] = np.round(np.median(sort_dratio_pool, axis = 0), d)
df["mean_density_ratio"] = np.round(np.mean(sort_dratio_pool, axis = 0), d)
df["std_density_ratio"] = np.round(np.std(sort_dratio_pool, axis = 0), d)
df["est_std_density_ratio"] = np.round(mad(sort_dratio_pool, axis = 0)/0.6745, d)
df["median_count_ratio"] = np.round(np.median(sort_cratio_pool, axis = 0), d)
df["mean_count_ratio"] = np.round(np.mean(sort_cratio_pool, axis = 0), d)
df["std_count_ratio"] = np.round(np.std(sort_cratio_pool, axis = 0), d)
df["est_std_count_ratio"] = np.round(mad(sort_cratio_pool, axis = 0)/0.6745, d)

df.index = grps

df.to_csv(os.path.join(sv_dst, "thal_ratio_stats.csv"))
