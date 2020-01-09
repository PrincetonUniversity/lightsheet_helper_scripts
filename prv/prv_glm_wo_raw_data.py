#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:02:15 2019

@author: wanglab
"""

import matplotlib as mpl, os
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
#imports
#path to pickle file
data_pth = "/jukebox/wang/zahra/tracing_projects/prv/for_tp/model_data_contra_pma.p"
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#set dest
dst = "/home/wanglab/Desktop"
#dst = "/Users/tjp7rr1/Downloads"

#set the appropritate variables
c_mat = data["c_mat"]
mat = data["mat"]
pmat = np.asarray(data["pmat"])
p_shuf = np.asarray(data["p_shuf"])
mat_shuf = np.asarray(data["mat_shuf"])
ak_pool = data["ak_pool"]
regions = data["regions"]
primary_lob_n = data["primary_lob_n"]

## display
fig = plt.figure(figsize=(6,5))
ax = fig.add_axes([.4,.1,.5,.8])

#set white text limit here
whitetext = 8
annotation_size = "medium"#annotation/number sizes

# map 1: weights
show = np.flipud(mat) # NOTE abs

vmin = 0
vmax = 10
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,((vmax-vmin)/2)+1)
#bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, 
                  boundaries=bounds, format="%d", shrink=0.3, aspect=10)
cb.set_label("Weight / SE", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

#annotations
for ri,row in enumerate(show): 
    for ci,col in enumerate(row):
        if col > whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize=annotation_size)
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize=annotation_size)   
            
# signif
sig = np.flipud(pmat) < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in np.flipud(regions)], fontsize="medium")

plt.savefig(os.path.join(dst, "prv_nc_glm_contra_layer56_pma.pdf"), bbox_inches = "tight")