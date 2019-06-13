#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 18:49:42 2019

@author: wanglab
"""

import numpy as np, matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#looking at goodness of fit for different models
ars_pc = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/r2adjusted_pc_layer.npy")
ars_pc_pool = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/r2adjusted_pc_layer_pooled_inj.npy")
ars = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/r2adjusted.npy")
ars_pool = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/r2adjusted_pooled_inj.npy")
regions = np.load("/jukebox/wang/zahra/modeling/h129/neocortex/neocortical_regions.npy")


fig, ax = plt.subplots(figsize=(13,10))

ax.scatter(range(17), ars_pc, color = "blue")
ax.axhline(0, color = 'grey')
ax.set_yticks(np.arange(-.5, .5, 10))
ax.set_xticks(np.arange(len(regions)))
ax.set_xticklabels(regions, rotation=20, fontsize=8, ha="right")
ax.set_xlabel("Neocortical regions")
ax.set_ylabel("Adjusted R-squared")

ax.scatter(range(17), ars_pc_pool, color = "green")
ax.axhline(0, color = 'grey')
ax.set_xticks(np.arange(len(regions)))
ax.set_xticklabels(regions, rotation=20, fontsize=8, ha="right")
ax.set_xlabel("Neocortical regions")
ax.set_ylabel("Adjusted R-squared")

ax.scatter(range(17), ars, color = "red")
ax.axhline(0, color = 'grey')
ax.set_xticks(np.arange(len(regions)))
ax.set_xticklabels(regions, rotation=20, fontsize=8, ha="right")
ax.set_xlabel("Neocortical regions")
ax.set_ylabel("Adjusted R-squared")

ax.scatter(range(17), ars_pool, color = "gold")
ax.axhline(0, color = 'grey')
ax.set_xticks(np.arange(len(regions)))
ax.set_xticklabels(regions, rotation=20, fontsize=8, ha="right")
ax.set_xlabel("Neocortical regions")
ax.set_ylabel("Adjusted R-squared")

ax.set_yticks(np.arange(-0.5, 0.6, 0.1))
ax.set_xticks(np.arange(len(regions)))

green_patch = mpatches.Patch(color="green", label="PC layer, pooled injections")
gold_patch = mpatches.Patch(color="gold", label="Entire regions, pooled injections")
blue_patch = mpatches.Patch(color="blue", label="PC layer")
red_patch = mpatches.Patch(color="red", label="Entire regions")

plt.legend(handles=[green_patch, gold_patch, blue_patch, red_patch], bbox_to_anchor=(.8, 1), loc=2, borderaxespad=0.)

plt.savefig('/home/wanglab/Desktop/r2sdj.svg')