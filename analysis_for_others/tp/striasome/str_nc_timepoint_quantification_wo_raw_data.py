#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 19:37:02 2019

@author: wanglab
"""


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, pandas as pd, seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#import data
pth = "/jukebox/wang/zahra/modeling/h129/striatum/count_and_density_data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

#set dst 
dst = "/home/wanglab/Desktop"
#dst = '/Users/tjp7rr1/Downloads/'

cell_counts_per_brain_p = data["cell_counts_per_brain_p"]
sois = data["sois"]
brains = data["brainnames"]
ak_pool = data["ak_pool"]
primary_pool = data["primary_pool"]
density_per_brain = data["density_per_brain"]
primary_lob_n = data["primary_lob_n"]

#mean percent counts
mean_counts = np.asarray([np.mean(cell_counts_per_brain_p[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,4))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T 

vmin = 0
vmax = 20
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
bounds = np.linspace(vmin,vmax,5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, shrink=0.5, aspect=10, format="%d")
cb.set_label("% of striatum counts", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
        
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="medium")
plt.savefig(os.path.join(dst,"str_mean_percent_counts.pdf"), bbox_inches = "tight")

mean_counts = np.asarray([np.mean(density_per_brain[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,4))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T 

vmin = 0
vmax = 300
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, shrink=0.5, aspect=10, format="%d")
cb.set_label("Cells/$mm^3$", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in sois], fontsize="medium")
plt.savefig(os.path.join(dst,"str_mean_density.pdf"), bbox_inches = "tight")

#%%
#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
order = np.argsort(np.mean(density_per_brain, axis = 0))[::-1]
sois_sort_density = np.array(sois)[order]

order = np.argsort(np.mean(cell_counts_per_brain_p, axis = 0))[::-1]
sois_sort_pcounts = np.array(sois)[order]

#boxplots of densitys
plt.figure()
df = pd.DataFrame(density_per_brain)
df.columns = sois 
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort_density)
sns.boxplot(data = df, orient = "h", showfliers=False,showcaps=False, boxprops={'facecolor':'None'}, order = sois_sort_density)
plt.xlabel("Density (cells/$mm^3$)")
plt.ylabel("Striatum structures")
plt.savefig(os.path.join(dst, "str_density_boxplots.pdf"), bbox_inches = "tight")

#boxplots of percent counts
plt.figure()
df = pd.DataFrame(cell_counts_per_brain_p)
df.columns = sois 
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort_pcounts)
sns.boxplot(data = df, orient = "h", showfliers=False,showcaps=False, boxprops={'facecolor':'None'}, order = sois_sort_pcounts)
plt.xlabel("% of total striatum cells")
plt.ylabel("Striatum structures")
plt.savefig(os.path.join(dst, "str_pcounts_boxplots.pdf"), bbox_inches = "tight")