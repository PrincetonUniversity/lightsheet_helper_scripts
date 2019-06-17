#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:02:47 2019

@author: wanglab
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt, matplotlib as mpl, os

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#set destination for figs
dst = "/home/wanglab/Desktop"

#import
pth = "/home/wanglab/Desktop/select_structures_percent_counts_for_visualization.csv"
df = pd.read_csv(pth)

#sort by animal
df = df.sort_values(by = ["animal"])

#remove homecage control for now
df = df[df.group != "homecage_control"]

#sort 1 - remove extremes
counts = np.asarray([df[xx].values for xx in df.columns if xx != "animal" and xx != "group" and np.mean(df[xx].values) < 3
                     and np.mean(df[xx].values) > 5e-1])

regions = np.asarray([xx for xx in df.columns if xx != "animal" and xx != "group" and np.mean(df[xx].values) < 3
                     and np.mean(df[xx].values) > 5e-1])

brains = df.animal.values

#formatting
fig = plt.figure(figsize=(11,4))
ax = fig.add_axes([.4,.1,.5,.7])

show = counts

vmin = 0
vmax = 3#.6
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,3,7)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("% of total counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="xx-small")
plt.savefig(os.path.join(dst,"cfos_percent_counts_1.pdf"), bbox_inches = "tight")

#%%

#sort 2 - high count regions
counts = np.asarray([df[xx].values for xx in df.columns if xx != "animal" and xx != "group" 
                     and np.mean(df[xx].values) > 3])

regions = np.asarray([xx for xx in df.columns if xx != "animal" and xx != "group" 
                     and np.mean(df[xx].values) > 3])
#formatting
fig = plt.figure(figsize=(11,1.2))
ax = fig.add_axes([.4,.1,.5,.7])

show = counts

vmin = 2
vmax = 10#.6
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,10,9)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.9, aspect=10)
cb.set_label("% of total counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="xx-small")

plt.savefig(os.path.join(dst,"cfos_percent_counts_2.pdf"), bbox_inches = "tight")

#%%

#sort 3 - low count regions
counts = np.asarray([df[xx].values for xx in df.columns if xx != "animal" and xx != "group" 
                     and np.mean(df[xx].values) < 5e-1])

regions = np.asarray([xx for xx in df.columns if xx != "animal" and xx != "group" 
                     and np.mean(df[xx].values) < 5e-1 ])
#formatting
fig = plt.figure(figsize=(11,4))
ax = fig.add_axes([.4,.1,.5,.7])

show = counts

vmin = 0
vmax = 0.5#.6
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("% of total counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="xx-small")

plt.savefig(os.path.join(dst,"cfos_percent_counts_3.pdf"), bbox_inches = "tight")

#%%

#sort 4 - dcn
counts = np.asarray([df[xx].values for xx in df.columns if xx == "Dentate nucleus" or xx == "Interposed nucleus"
                     or xx == "Fastigial nucleus"])

regions = np.asarray([xx for xx in df.columns if xx == "Dentate nucleus" or xx == "Interposed nucleus"
                     or xx == "Fastigial nucleus"])
#formatting
fig = plt.figure(figsize=(11,0.5))
ax = fig.add_axes([.4,.1,.5,.8])

show = counts

vmin = 0
vmax = 0.15#.6
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.3f", 
                  shrink=0.3, aspect=10)
cb.set_label("% of total counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="xx-small")

plt.savefig(os.path.join(dst,"cfos_percent_counts_dcn.pdf"), bbox_inches = "tight")

#%%
#get and plot behavioral metrics
behav = "/home/wanglab/Desktop/ymaze.csv"

bh = pd.read_csv(behav)

bhm = np.asarray([bh[xx].values for xx in bh.columns if xx == "Initial Reversal" or xx == "Final Reversal"
                  or xx == "Forced Reversal"])

bhn = np.asarray([xx for xx in bh.columns if xx == "Initial Reversal" or xx == "Final Reversal"
                  or xx == "Forced Reversal"])

#plot reversal percents
#formatting
fig = plt.figure(figsize=(11,.5))
ax = fig.add_axes([.4,.1,.5,.7])

show = bhm

vmin = 0
vmax = 100#.6
cmap = plt.cm.Reds
cmap.set_over('maroon')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=1, aspect=10)
cb.set_label("Performance", fontsize=3, labelpad=3)
cb.ax.tick_params(labelsize=3)

cb.ax.set_visible(True)
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(bhn))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in bhn], fontsize="xx-small")
plt.savefig(os.path.join(dst,"behav.pdf"), bbox_inches = "tight")

#%%
#plot multisession reversal + learning
bhm = np.asarray([bh[xx].values for xx in bh.columns if xx == "Multisession Reversal " or xx == "Multisession Learning"])

bhn = np.asarray([xx for xx in bh.columns if xx == "Multisession Reversal " or xx == "Multisession Learning"])

fig = plt.figure(figsize=(11,.4))
ax = fig.add_axes([.4,.1,.5,.7])

show = bhm

vmin = -0.02
vmax = 0.1#.6
cmap = plt.cm.Reds
cmap.set_over('maroon')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.3f", 
                  shrink=1, aspect=10)
cb.set_label("Performance", fontsize=3, labelpad=3)
cb.ax.tick_params(labelsize=3)

cb.ax.set_visible(True)
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(bhn))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in bhn], fontsize="xx-small")
plt.savefig(os.path.join(dst,"behav_multisession_reversal.pdf"), bbox_inches = "tight")

#%%
#plot distance traveled
bhm = np.asarray([bh["Ymaze Distance Traveled"].values, bh["Ymaze Distance Traveled"].values])

bhn = np.asarray(["Ymaze Distance Traveled", "Ymaze Distance Traveled"])

fig = plt.figure(figsize=(11,.2))
ax = fig.add_axes([.4,.1,.5,.7])

show = bhm

vmin = 8500
vmax = 15000#.6
cmap = plt.cm.Reds
cmap.set_over('maroon')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=3, aspect=10)
cb.set_label("Distance", fontsize=3, labelpad=3)
cb.ax.tick_params(labelsize=3)

cb.ax.set_visible(True)
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(1)+1)
ax.set_yticklabels(["Ymaze Distance Traveled"], fontsize="xx-small")
plt.savefig(os.path.join(dst,"behav_distance.pdf"), bbox_inches = "tight")

#%%

#plot inj fractions

#injection columns
inj = ['DREADD voxel fraction', 'lobvi', 'lobvii', 'crus1', 'crus2', 'simplex', 'vermis',
       'hemisphere']

bhm = np.asarray([bh[xx].values for xx in bh.columns if xx in inj])

bhn = np.asarray(['cerebellar fraction', 'Lob. VI', 'Lob. VII', 'Crus 1', 'Crus 2', 'Simplex',
                  'Vermis', 'Hemisphere', ])


fig = plt.figure(figsize=(11,2))
ax = fig.add_axes([.4,.1,.5,.7])

show = bhm

vmin = 0
vmax = 0.4#.6
cmap = plt.cm.Blues
cmap.set_over('lightslategray')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,5)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=.3, aspect=10)
cb.set_label("Distance", fontsize=3, labelpad=3)
cb.ax.tick_params(labelsize=3)

cb.ax.set_visible(True)

ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(bhn))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in bhn], fontsize="xx-small")

plt.savefig(os.path.join(dst,"inj.pdf"), bbox_inches = "tight")