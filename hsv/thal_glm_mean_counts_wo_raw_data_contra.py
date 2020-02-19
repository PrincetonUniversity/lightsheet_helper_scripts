#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 17:50:18 2019

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#imports
#path to pickle file
data_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_model_data_contra_allen.p"

data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#set dest
dst = "/home/wanglab/Desktop"

#set the appropritate variables
c_mat = data["c_mat"]
mat = data["mat"]
pmat = np.asarray(data["pmat"])
p_shuf = np.asarray(data["p_shuf"])
mat_shuf = np.asarray(data["mat_shuf"])
ak_pool = data["ak_pool"]
regions = data["regions"]
primary_lob_n = data["primary_lob_n"]
primary_pool = data["primary_pool"]
#get cell counts from raw data export?
pcounts = data["pcounts"]


#change nuclei names (for y axis)here!!
regions = np.array(["VPM", "VPL","VA-L", "AV",  "LD",
       "PV", "MedHab", "LP",
       "PoT", "MD", "Po",
       "VM", "RTN"])

#glm figure
## display
fig = plt.figure(figsize=(5.5,5))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat

vmin = 0
vmax = 5
whitetext = 4
cmap = plt.cm.Blues
cmap.set_under("w")
cmap.set_over(plt.cm.Blues(1.0))#cmap.set_over("maroon")
annotation = False
#colormap

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Model weight / SE", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
if annotation:
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
    ax.text(x, y+0.2, "*", fontsize=12, ha="left", va="bottom", color = "white", transform=ax.transData)

# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(ak_pool)
ax.set_xticklabels([])#["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels([])#["{}".format(bi) for bi in regions], fontsize="small")

#despline to make it look similar to paper figure
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.grid(False)
ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")

#%%

# SET COLORMAP
cmap = plt.cm.Greys
cmap.set_over(cmap(1.0))

#set min and max of colorbar
vmin = 0
vmax = 8

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(pcounts[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5.5,5))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Mean % of thalamic counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels([])#["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels([])#["{}".format(bi) for bi in regions], fontsize="small")

#despline to make it look similar to paper figure
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.grid(False)
ax.tick_params(length=6)

plt.savefig(os.path.join(dst,"thal_mean_count.pdf"), bbox_inches = "tight")
