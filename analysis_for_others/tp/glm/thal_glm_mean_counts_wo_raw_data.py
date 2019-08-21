#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 13:38:57 2019

@author: wanglab
"""

import matplotlib as mpl, os
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

#imports
#path to pickle file
data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/model_data_v2.p"
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#raw data
raw_data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data_v2.p"
raw = pckl.load(open(raw_data_pth, "rb"), encoding = "latin1")

#set dest
dst = "/home/wanglab/Desktop"
#dst = "/Users/tjp7rr1/Downloads"

#set the appropritate variables
c_mat = data["c_mat"]
mat = data["mat"]
pmat = np.asarray(data["pmat"])
p_shuf = np.asarray(data["p_shuf"])
mat_shuf = np.asarray(data["mat_shuf"])
ak = data["ak_pool"]
regions = data["regions"]
primary_lob_n = data["primary_lob_n"]
primary_pool = raw["primary_pool"]
#get cell counts from raw data export?
cell_counts_per_brain_p = raw["cell_counts_per_brain"]

#set to reorder based on jones classification
re_ord = True

if re_ord:
    new_ord = np.array(["Lateral habenula", 
                          "Lateral posterior nucleus of the thalamus",
                          "Lateral dorsal nucleus of thalamus",
                          "Posterior triangular thalamic nucleus",
                          "Posterior complex of the thalamus",
                          "Parafascicular nucleus", 
                          "Central lateral nucleus of the thalamus",
                          "Paraventricular nucleus of the thalamus", "Nucleus of reuniens",
                          "Submedial nucleus of the thalamus",
       "Mediodorsal nucleus of thalamus",
       "Ventral part of the lateral geniculate complex",
       "Reticular nucleus of the thalamus",
       "Anteroventral nucleus of thalamus",
       "Ventral medial nucleus of the thalamus",
       "Ventral anterior-lateral complex of the thalamus",
       "Ventral posterolateral nucleus of the thalamus",
       "Ventral posteromedial nucleus of the thalamus"])
    ord_idx = [i for j, yy in enumerate(new_ord) for i, xx in enumerate(regions) if new_ord[j] == regions[i]]
    regions = regions[ord_idx]
    mat = mat[ord_idx]
    pmat = pmat[ord_idx]
    cell_counts_per_brain_p = cell_counts_per_brain_p.T[ord_idx]

#show actual numbers for % counts for the mean counts figure
show_ann = False

#glm figure
## display
fig = plt.figure(figsize=(5,7))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
#changed order of matrix, so that regions are reversed
show = np.flipud(mat)#np.median(mat_shuf, axis = 0) # NOTE abs

vmin = 0
vmax = 3
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.2, aspect=10)
cb.set_label("Weight / SE", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(c_mat):
    for ci,col in enumerate(row):
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# signif
sig = pmat < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
#flipped matrix with signifcant pvals too 
for y,x in np.argwhere(np.flipud(sig)):
    pass
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05", ha="center", va="center", fontsize="small", transform=ax.transAxes)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *"s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
# xticksjt -t monokai -m 200
ax.set_xticks(np.arange(len(ak))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(ak)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
#changing order of regions for tp to visualize
ax.set_yticklabels(["{}".format(bi) for bi in np.flip(regions)], fontsize="xx-small")
plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")


#mean counts
#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(cell_counts_per_brain_p[:,np.where(primary_pool == idx)[0]], axis=1) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,7))
ax = fig.add_axes([.4,.1,.5,.8])

show = np.flipud(mean_counts.T) 

vmin = 0
vmax = 6
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("Mean % of thalamic counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
if show_ann:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            pass
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")
        
# aesthetics
ax.set_xticks(np.arange(len(ak))+.5)
lbls = np.asarray(ak)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
#changing order of regions for tp to visualize
ax.set_yticklabels(["{}".format(bi) for bi in np.flip(regions)], fontsize="xx-small")
plt.savefig(os.path.join(dst, "thal_mean_counts.pdf"), bbox_inches = "tight")