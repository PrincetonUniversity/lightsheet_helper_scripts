#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 19:37:02 2019

@author: wanglab
"""


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#import data
pth = "/jukebox/wang/zahra/modeling/h129/striatum/count_and_density_data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

#set dst 
dst = "/home/wanglab/Desktop"

cell_counts_per_brain_p = data["cell_counts_per_brain_p"]
sois = data["sois"]
brains = data["brainnames"]
ak_pool = data["ak_pool"]
primary_pool = data["primary_pool"]
density_per_brain = data["density_per_brain"]
primary_lob_n = data["primary_lob_n"]

#mean percent counts
mean_counts = np.asarray([np.mean(cell_counts_per_brain_p[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(4,3))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T 

vmin = 0
vmax = 20
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,5)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.3, aspect=10)
cb.set_label("% of striatum counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 3:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="xx-small")
plt.savefig(os.path.join(dst,"str_mean_percent_counts.pdf"), bbox_inches = "tight")

mean_counts = np.asarray([np.mean(density_per_brain[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(4,3))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T 

vmin = 0
vmax = 300
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,4)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.3, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 50:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="xx-small")
plt.savefig(os.path.join(dst,"str_mean_density.pdf"), bbox_inches = "tight")


#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
density_per_brain_descending_order = np.sort(density_per_brain)
order = np.argsort(np.mean(density_per_brain, axis = 0))
sois_descending_density = np.array(sois)[order]

cell_counts_per_brain_p_descending_order = np.sort(cell_counts_per_brain_p)
order = np.argsort(np.mean(cell_counts_per_brain_p, axis = 0))
sois_descending_pcounts = np.array(sois)[order]

#boxplots of densitys
plt.figure(figsize=(5,6))
plt.boxplot(density_per_brain_descending_order, vert = False, labels = sois_descending_density, sym = "", showcaps = False)
ngroup = len(density_per_brain_descending_order.T)
for i in range(ngroup):
    plt.scatter(density_per_brain_descending_order[:,i], y=np.ones(len(density_per_brain_descending_order[:,0]))*i+1, color = "k", s = 7)
plt.xlabel("Density (cells/$mm^3$)")
plt.ylabel("Striatum structures")
plt.savefig(os.path.join(dst, "str_density_boxplots.pdf"), bbox_inches = "tight")

#boxplots of percent counts
plt.figure(figsize=(5,6))
plt.boxplot(cell_counts_per_brain_p_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(cell_counts_per_brain_p_descending_order.T)
for i in range(ngroup):
    plt.scatter(cell_counts_per_brain_p_descending_order[:,i], 
                y=np.ones(len(cell_counts_per_brain_p_descending_order[:,0]))*i+1, color = "k", s = 10)
plt.xlabel("% of total striatum cells")
plt.ylabel("Striatum structures")
plt.savefig(os.path.join(dst, "str_pcounts_boxplots.pdf"), bbox_inches = "tight")