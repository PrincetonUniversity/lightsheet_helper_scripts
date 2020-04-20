#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:44:22 2019

@author: wanglab
"""

import matplotlib as mpl, itertools
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#import data
pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/vtasnc_counts_contra.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

#set dst 
dst = "/home/wanglab/Desktop"

counts_per_struct = data["counts_per_struct"]
frac_of_inj_pool = data["frac_of_inj_pool"]
brains = data["brainnames"]
ak_pool = data["ak_pool"]
primary_pool = data["primary_pool"]
short_nuclei = data["short_nuclei"]
density_per_struct = data["density_per_struct"]
primary_lob_n = data["primary_lob_n"]

#%%

#INJ FRACTION MAP

fig, ax = plt.subplots(figsize = (10,2))

sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
show = np.fliplr(sort_inj).T

vmin = 0.05
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
cmap.set_under('white')

#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("% coverage of lobule", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize="xx-small", ha="right")

#despline to make it look similar to paper figure
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.grid(False)

plt.savefig(os.path.join(dst, "vta_inj.pdf"), bbox_inches = "tight")

#%%
## CELL COUNTS
   
#sort 
sorted_counts = [counts_per_struct.T[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sorted_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]

#reformat - wtf
import itertools
sorted_counts = np.asarray(list(itertools.chain.from_iterable(sorted_counts)))
sorted_brains = list(itertools.chain.from_iterable(sorted_brains))
    
fig = plt.figure(figsize=(15,2))
ax = fig.add_axes([.4,.1,.5,.8])

show = sorted_counts.T

vmin = 0
vmax = 75
whitetext = 20
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.8, aspect=12)
cb.set_label("Cell counts", fontsize="small", labelpad=2)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="medium")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="medium")
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(sorted_brains, rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="medium")

#despline to make it look similar to paper figure
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.grid(False)

plt.savefig(os.path.join(dst,"thalvtacomp_cell_counts_sorted.pdf"), bbox_inches = "tight")

#%%
## DENSITY
#sort density
sorted_counts = [density_per_struct[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sorted_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
import itertools
sorted_counts = np.asarray(list(itertools.chain.from_iterable(sorted_counts)))
sorted_brains = list(itertools.chain.from_iterable(sorted_brains))
    
fig = plt.figure(figsize=(15,2))
ax = fig.add_axes([.4,.1,.5,.8])

show = sorted_counts.T.astype(int) #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 75
whitetext = 20
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.8, aspect=12)
cb.set_label("Cells/$mm^3$", fontsize="small", labelpad=2)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="medium")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="medium")
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(sorted_brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(sorted_brains, rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="medium")

#despline to make it look similar to paper figure
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.grid(False)

plt.savefig(os.path.join(dst,"thalvtacomp_density_sorted.pdf"), bbox_inches = "tight")
    
#%%
## mean density
fig = plt.figure(figsize=(5,2))
ax = fig.add_axes([.4,.1,.5,.8])

mean_counts = np.asarray([np.mean(density_per_struct[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

show = mean_counts.T.astype(int) #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 75
whitetext = 20
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.8, aspect=12)

cb.set_label("Cells/$mm^3$", fontsize="small", labelpad=2)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="medium")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="medium")
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(ak_pool, primary_lob_n)], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="medium")

#despline to make it look similar to paper figure
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.grid(False)

plt.savefig(os.path.join(dst,"thalvtacomp_mean_density.pdf"), bbox_inches = "tight")
#%%

## mean counts
fig = plt.figure(figsize=(5,2))
ax = fig.add_axes([.4,.1,.5,.8])

mean_counts = np.asarray([np.mean(counts_per_struct.T[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

show = mean_counts.astype(int).T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 75
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.8, aspect=12)

cb.set_label("Count", fontsize="small", labelpad=2)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < 30:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="medium")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="medium")
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(ak_pool, primary_lob_n)], rotation=45, fontsize=7, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="medium")

#despline to make it look similar to paper figure
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.grid(False)

plt.savefig(os.path.join(dst,"thalvtacomp_mean_counts.pdf"), bbox_inches = "tight")