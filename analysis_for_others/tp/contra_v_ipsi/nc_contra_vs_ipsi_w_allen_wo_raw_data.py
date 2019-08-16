#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 19:23:50 2019

@author: wanglab
"""


import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, matplotlib as mpl, pandas as pd
from scipy.stats import median_absolute_deviation as mad

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#import data
main_data_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_contra_ipsi_counts_densities.p"
data = pckl.load(open(main_data_pth, "rb"), encoding = "latin1")

#injection site analysis
inj_pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
inj_dct = pckl.load(open(inj_pth, "rb"), encoding = "latin1")

#set save destination
sv_dst = "/home/wanglab/Desktop"

cell_counts_per_brain_left = data["cell_counts_per_brain_left"]
cell_counts_per_brain_right = data["cell_counts_per_brain_right"]
density_per_brain_left = data["density_per_brain_left"]
density_per_brain_right = data["density_per_brain_right"]
lr_dist = data["lr_dist"]
nc_areas = data["nc_areas"] #gives order of nc areas also

brains = inj_dct["brainnames"]
primary_pool = inj_dct["primary_pool"]
ak_pool = inj_dct["cb_regions_pool"]
inj = inj_dct["expr_all_as_frac_of_inj_pool"]

#-------------------------------------------------------------------------------------------------------------------------------------
#preprocessing
nc_left_counts = cell_counts_per_brain_left
nc_right_counts = cell_counts_per_brain_right
nc_density_left = density_per_brain_left
nc_density_right = density_per_brain_right

lrv = list(lr_dist.values())
lr_brains = list(lr_dist.keys())

#dct is just for my sanity, so im not mixing up brains
_ccontra = []; _cipsi = []; _dcontra = []; _dipsi = []
for i in range(len(lr_brains)):
    if lrv[i] > 0: #right
        #counts
        _ccontra.append(nc_left_counts[i])
        _cipsi.append(nc_right_counts[i])
        #density
        _dcontra.append(nc_density_left[i])
        _dipsi.append(nc_density_right[i])
    elif lrv[i] < 0: #left
        #counts
        _ccontra.append(nc_right_counts[i])
        _cipsi.append(nc_left_counts[i])
        #density
        _dcontra.append(nc_density_right[i])
        _dipsi.append(nc_density_left[i])


_ccontra = np.asarray(_ccontra).T; _dcontra = np.asarray(_dcontra).T
_cipsi = np.asarray(_cipsi).T; _dipsi = np.asarray(_dipsi).T
_dratio = np.asarray([_dcontra[i]/_dipsi[i] for i in range(len(_dcontra))])
_cratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
#make into one
_dist = np.asarray(list(lr_dist.values()))

#injection site analysis
pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

brains = data["brainnames"]
primary_pool = data["primary_pool"]
ak_pool = data["cb_regions_pool"]
inj = data["expr_all_as_frac_of_inj_pool"]
 
_inj = np.asarray([inj[i] for i in range(len(inj)) if brains[i] in lr_brains])
_primary_pool = np.asarray([primary_pool[i] for i in range(len(primary_pool)) if brains[i] in lr_brains])

#sort by distance
sort_dist = np.sort(_dist)
sort_ccontra = _ccontra.T[np.argsort(_dist, axis = 0)]
sort_cipsi = _cipsi.T[np.argsort(_dist, axis = 0)]
sort_cratio = _cratio.T[np.argsort(_dist, axis = 0)]
sort_dcontra = _dcontra.T[np.argsort(_dist, axis = 0)]
sort_dipsi = _dipsi.T[np.argsort(_dist, axis = 0)]
sort_dratio = _dratio.T[np.argsort(_dist, axis = 0)]

sort_inj = _inj[np.argsort(_dist)]   
sort_brains = np.array(lr_brains)[np.argsort(_dist)]

print(sort_dist.shape)
print(sort_cratio.shape)

#group into 3 regions
grps = np.array(["Frontal\n(IL,PL,ACC,ORB,FRP,\nGU,AI,VISC)" , "Medial\n(MO,SS)", "Posterior\n(RSP,PTL,TE,PERI,ECT)"])
sort_ccontra_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_ccontra])
sort_dcontra_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_dcontra])
sort_cipsi_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_cipsi])
sort_dipsi_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_dipsi])
sort_cratio_pool = np.asarray([sort_ccontra_pool[i]/sort_cipsi_pool[i] for i in range(len(sort_ccontra_pool))])
sort_dratio_pool = np.asarray([sort_dcontra_pool[i]/sort_dipsi_pool[i] for i in range(len(sort_dcontra_pool))])

#-------------------------------------------------------------------------------------------------------------------------------------

## display
fig, axes = plt.subplots(ncols = 1, nrows = 5, figsize = (15,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,1,1,1,0.3]})

#set colormap specs
vmaxcount = 500
whitetext = 200
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

ax = axes[1]
show = sort_dcontra_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcount
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
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
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[2]
show = sort_dipsi_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcount
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
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

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)


ax = axes[3]
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

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra/Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[4]
show = np.asarray([sort_dist])

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

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

plt.savefig("/home/wanglab/Desktop/nc_density_ratios.pdf", bbox_inches = "tight")

#-------------------------------------------------------------------------------------------------------------------------------------
## display
fig, axes = plt.subplots(ncols = 1, nrows = 5, figsize = (15,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,1,1,1,0.3]})

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

ax = axes[1]
show = sort_ccontra_pool.T
yaxis = grps

vmin = 0
vmax = 2000
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
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 200:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[2]
show = sort_cipsi_pool.T
yaxis = grps

vmin = 0
vmax = 2000
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
        if col < 200:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)


ax = axes[3]
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
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[4]
show = np.asarray([sort_dist])

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
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=5, ha="right")

ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

plt.savefig("/home/wanglab/Desktop/nc_count_ratios.pdf", bbox_inches = "tight")

#basic stats for these ratios
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

df.to_csv(os.path.join(sv_dst, "nc_ratio_stats.csv"))