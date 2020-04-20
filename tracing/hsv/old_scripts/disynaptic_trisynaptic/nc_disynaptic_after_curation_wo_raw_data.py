#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:22:57 2019

@author: wanglab
"""

import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, matplotlib as mpl, statsmodels.api as sm, pandas as pd
from scipy.stats import median_absolute_deviation as mad

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#import data
pth = "/jukebox/wang/zahra/h129_qc/data/nc_density_at_thal_nc_timepoint_data_all_brains.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

#set dst 
dst = "/home/wanglab/Desktop"

#mask to remove sandy brains
curated_brains = [False, True, True, False, False, False, True, False, True, False, True, True, False, False, 
                  False, False, True, False, True, False, True, False, True]

thal_density_per_brain = data["thal_density_per_brain"][curated_brains]
lbls = data["nc_lbls"]
nc_brains = data["nc_brainnames"]
thal_brains = data["thal_brainnames"][curated_brains]
nc_density_per_brain = data["nc_density_per_brain"]
mean_thal_density_per_brain = np.mean(thal_density_per_brain, axis=0)
mean_nc_density_per_brain = np.mean(nc_density_per_brain, axis=0)
std_thal_density_per_brain = np.std(thal_density_per_brain, axis=0)
std_nc_density_per_brain = np.std(nc_density_per_brain, axis=0)

#-----------------------------------------------------------------------------------------------------------------------------
## display
fig = plt.figure(figsize=(10, 5))
ax = fig.add_axes([.4,.1,.5,.8])

show = np.flipud(thal_density_per_brain.T)

vmin = 0
vmax = 100
cmap = plt.cm.viridis
cmap.set_under("w")
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.2, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 20:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# aesthetics
ax.set_xticks(np.arange(len(thal_brains))+.5)
ax.set_xticklabels(np.asarray(thal_brains), rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(lbls))+.5)
ax.set_yticklabels(np.flipud(np.asarray(lbls)), fontsize="x-small")
ax.set_ylabel("Thalamic 'disynaptic' timepoint", fontsize="x-small")
ax.yaxis.set_label_coords(-0.3,0.5)

plt.savefig(os.path.join(dst,"nc_density_at_thalamic_timepoint.pdf"), bbox_inches = "tight")

#-----------------------------------------------------------------------------------------------------------------------------
#NC
fig = plt.figure(figsize=(22, 5))
ax = fig.add_axes([.4,.1,.5,.8])

show = np.flipud(nc_density_per_brain.T)

vmin = 0
vmax = 500
cmap = plt.cm.viridis
cmap.set_under("w")
cmap.set_over("gold")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.2, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 100:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(nc_brains))+.5)
#remaking labeles so it doesn"t look squished
nc_brains = np.asarray(nc_brains)
ax.set_xticklabels(nc_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(lbls))+.5)
ax.set_yticklabels(np.flipud(np.asarray(lbls)), fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Neocortical 'trisynaptic' timepoint", fontsize="x-small")
ax.yaxis.set_label_coords(-0.22,0.5)

plt.savefig(os.path.join(dst,"nc_density_at_nc_timepoint.pdf"), bbox_inches = "tight")

#-----------------------------------------------------------------------------------------------------------------------------
#save the stats for the densities into csv
est_std = 0.6745 #normal mad to get estimated standard dev
ratio_mean_density = np.array(mean_thal_density_per_brain/mean_nc_density_per_brain)
ratio_std_density = np.array(std_thal_density_per_brain/std_nc_density_per_brain)
#calculate median also
median_thal_density_per_brain = np.median(thal_density_per_brain, axis = 0)
median_nc_density_per_brain = np.median(nc_density_per_brain, axis = 0)
ratio_median_density = np.array(median_thal_density_per_brain/median_nc_density_per_brain)

mad_thal_density_per_brain = mad(thal_density_per_brain, axis = 0)
mad_nc_density_per_brain = mad(nc_density_per_brain, axis = 0)
ratio_mad_density = np.array(mad_thal_density_per_brain/mad_nc_density_per_brain)

df = pd.DataFrame()
d = 4 #decimals to round to
df["mean_thal_density"] = np.round(mean_thal_density_per_brain, d)
df["mean_nc_density"] = np.round(mean_nc_density_per_brain, d)
df["std_thal_density"] = np.round(std_thal_density_per_brain, d)
df["std_nc_density"] = np.round(std_nc_density_per_brain, d)
df["median_thal_density"] = np.round(median_thal_density_per_brain, d)
df["median_nc_density"] = np.round(median_nc_density_per_brain, d)
df["est_std_thal_density"] = np.round(mad_thal_density_per_brain/est_std, d)
df["est_std_nc_density"] = np.round(mad_nc_density_per_brain/est_std, d)

df["ratio_mean_density"] = np.round(ratio_mean_density, d)
df["ratio_median_density"] = np.round(ratio_median_density, d)
df["ratio_std_density"] = np.round(ratio_std_density, d)
df["ratio_est_std_density"] = np.round(ratio_mad_density/est_std, d)
df.index = lbls

#make another dataframe with just ratio stats
df_stats = pd.DataFrame()
df_stats["mean"] = np.round(df.mean(axis = 0).values, d)
df_stats["median"] = np.round(df.median(axis = 0).values, d)
df_stats["standard_deviation"] = np.round(df.std(axis = 0).values, d)
df_stats["est_standard_deviation"] = np.round(mad(df.to_numpy(), axis = 0)/est_std, d)
df_stats.index = df.columns
df = df.append(df_stats.T)

df.to_csv(os.path.join(dst, "disynaptic.csv"))


#-----------------------------------------------------------------------------------------------------------------------------

#fig = plt.figure(figsize=(10,5))
#ax = fig.add_axes([.4,.1,.5,.8])
#
##linear regression - not useful
#Y = np.sort(mean_nc_density_per_brain)
#X = mean_thal_density_per_brain[np.argsort(mean_nc_density_per_brain)][:-1]
#strcs = np.asarray(lbls)[np.argsort(mean_nc_density_per_brain)][:-1]
#results = sm.OLS(Y,sm.add_constant(X)).fit()
#
#mean_slope = results.params[0]
#mean_r2 = results.rsquared
#mean_intercept = results.params[1]
##plot as scatter
##make them all different colors
#color_map = ["dimgray", "rosybrown", "darkred", "tomato", "chocolate", "orange", "gold", "olive", "darkseagreen", "springgreen", "teal",
#             "darkturquoise", "steelblue", "navy", "indigo", "crimson", "deeppink"]
#
#size = 30
#for i in range(len(X)):
#    ax.scatter(y = Y[i], x = X[i], s = size, label = strcs[i], c = color_map[i])
##
#X = range(0, 30)
#ax.plot(X, X*mean_slope + mean_intercept, "--r", label = "Slope={}\n$R^2$={}".format(round(mean_slope, 2), 
#                   round(mean_r2, 2)))
#
# 
#ax.set_ylim([0, 500])
#ax.set_xticks(np.arange(0, 30, 2))
##ax.set_xlim([0, 100])
#ax.set_ylabel("Average neocortical densities at neocortical timepoint")
#ax.set_xlabel("Average neocortical densities at thalamic timepoint")
#plt.legend(prop={'size': 10}, bbox_to_anchor=(1,1), loc='upper left', ncol=1)
#plt.savefig(os.path.join(dst, "disynaptic.pdf"), bbox_inches = "tight")