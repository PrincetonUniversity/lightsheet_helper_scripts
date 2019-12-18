#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 17:50:18 2019

@author: wanglab
"""

import matplotlib as mpl, os
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

#imports
#path to pickle file
data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/thal_model_data_contra_allen.p"
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
primary_pool = data["primary_pool"]
#get cell counts from raw data export?
pcounts = data["pcounts"]


#change nuclei names (for y axis)here!!
regions = np.array(['VPM', 'VPL','VA-L', 'AV',  'LD',
       'Paraventricular', 'M habenula', 'LP',
       'Post. Triangular', 'MD', 'Post. complex',
       'VM', 'RTN'])

#glm figure
## display
fig = plt.figure(figsize=(6,5))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat

vmin = 0
vmax = 5
whitetext = 4
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.2, aspect=10)
cb.set_label("Weight / SE", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
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
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05", ha="center", va="center", fontsize="small", transform=ax.transAxes)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="small")

plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")


#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(pcounts[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(6,5))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 8
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,(vmax-vmin)/2 +1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("Mean % of thalamic counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
## exact value annotations
#for ri,row in enumerate(show):
#    for ci,col in enumerate(row):
#        pass
#        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="small")
plt.savefig(os.path.join(dst,"thal_mean_count.pdf"), bbox_inches = "tight")
