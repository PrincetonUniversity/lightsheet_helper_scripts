#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 17:03:52 2020

@author: wanglab
"""


import matplotlib as mpl, os
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl


#TP
plt.rcParams["axes.grid"] = False

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
#imports
#path to pickle file
data_pth = "/home/wanglab/Desktop/zahra/nc_model_data_contra_pma.p"#/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_model_data_contra_pma.p"
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

vmin = 1
vmax = 7
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")


#tp local
tp = True
if tp:
    # dst = "/Users/tjp7rr1/Downloads"
    vmin = 0
    vmax = 5
    cmap = plt.cm.Blues#plt.cm.Reds
    cmap.set_over(plt.cm.Blues(1.0)) #cmap.set_over('maroon')
    whitetext = 4
    cmap.set_under('w')
    #reorder xaxis
    if False: #no longer changing the order - this will break too since new categories
        pmat = np.concatenate([pmat[:,:3], pmat[:,3:][:,::-1]], 1)
        ak_pool = np.concatenate([ak_pool[:3], ak_pool[3:][::-1]], 0)
        primary_lob_n = np.concatenate([primary_lob_n[:3], primary_lob_n[3:][::-1]], 0)
    
    #TP NOTE - FOR SOME REASON THIS MESSES UP SIGNIFICANT ASTERISKS, numbers are ok
    #LOOK AT ZAHRAS FOR THAT



#colormap
bounds = np.linspace(vmin,vmax,((vmax-vmin))+1)
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.3, aspect=10)
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

#despline to make it look similar to paper figure
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.grid(False)

plt.savefig(os.path.join(dst, "h129_nc_glm_contra_pma.pdf"), bbox_inches = "tight")
#%%