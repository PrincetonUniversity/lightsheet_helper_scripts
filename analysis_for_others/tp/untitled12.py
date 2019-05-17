#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:00:51 2019

@author: wanglab
"""

import pickle as pckl, numpy as np, seaborn as sns, pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline

#path to pickle file
data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data.p"
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")
print(data.keys())

#set the appropritate variables
X = data["expr_all_as_frac_of_lob_pool_norm"]
Y = data["cell_counts_per_brain"]
regions = data["thl_regions"] #for figures/annotations
ak = data["cb_regions_pool"]
#will use this for figure
primary_lob_n = np.asarray([np.where(data["primary_pool"] == i)[0].shape[0] for i in np.unique(data["primary_pool"])])
dst = "/home/wanglab/Desktop"

#%%
#organize this by injection?
primary = data["primary_pool"] 
Y_sort = np.asarray([Y[21], Y[0], Y[8], Y[14], Y[16], Y[18], Y[20], Y[11], Y[5], Y[13],
          Y[2], Y[3], Y[7], Y[10], Y[12], Y[15], Y[17], Y[19], Y[22],
          Y[1], Y[4], Y[6], Y[9]])

bn = data["brainnames"]
brainnames_sort = np.asarray([bn[21], bn[0], bn[8], bn[14], bn[16], bn[18], bn[20], bn[11], bn[5], bn[13],
          bn[2], bn[3], bn[7], bn[10], bn[12], bn[15], bn[17], bn[19], bn[22],
          bn[1], bn[4], bn[6], bn[9]])

fig, ax = plt.subplots(figsize = (8, 7))
plt.imshow(Y_sort.T)
ax.set_yticks(np.arange(len(regions)))
ax.set_yticklabels(["{}".format(ak) for ak in regions], fontsize=7, va = "center_baseline", ha = "right")
ax.set_xticks(np.arange(len(brainnames_sort)))
ax.set_xticklabels(["{}".format(brain) for brain in brainnames_sort], rotation=30, fontsize=5, ha="right")
cb = plt.colorbar(shrink = 0.3)
cb.set_label("% of thalamic counts", fontsize="x-small", labelpad=3)
plt.savefig(dst+"/counts_sort.svg", bbox_inches = "tight")

#%%
plt.style.use('default')
plt.figure(figsize = (12, 10))
#dendogram analysis
df = pd.DataFrame(Y, bn)
df.columns = ['Ventral group of the dorsal thalamus',
       'Subparafascicular nucleus', 'Subparafasicular area',
       'Peripeduncular nucleus', 'Geniculate group, dorsal thalamus',
       'Lateral group of the dorsal thalamus',
       'Anterior group of the dorsal thalamus',
       'Medial group of the dorsal thalamus',
       'Midline group of the dorsal thalamus',
       'Intralaminar nuclei of the dorsal thalamus',
       'Reticular nucleus of the thalamus',
       'Geniculate group, ventral thalamus', 'Epithalamus']

sns.clustermap(df.T, cmap = "viridis", metric="correlation")
sns.axes_style({ 'patch.force_edgecolor': False})
plt.savefig('/home/wanglab/Desktop/test.svg', bbox_inches = "tight")
#%%
X = data["expr_all_as_frac_of_lob_pool"]
X_sort = np.asarray([X[1], X[4], X[19], X[9], X[16], X[5],
                     X[22], X[8], X[11], X[18], X[12], 
                     X[2], X[7], X[6], X[21], X[15], X[20],
                     X[3], X[13], X[10], X[17], X[0], X[14]])

bn_sort = np.asarray([bn[1], bn[4], bn[19], bn[9], bn[16], bn[5],
                     bn[22], bn[8], bn[11], bn[18], bn[12], 
                     bn[2], bn[7], bn[6], bn[21], bn[15], bn[20],
                     bn[3], bn[13], bn[10], bn[17], bn[0], bn[14]])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "firebrick"]) 
plt.subplots(figsize = (6, 2.5))
inj = pd.DataFrame(X_sort, bn_sort)
inj.columns = ak
sns.heatmap(inj.T, cmap = cmap, cbar_kws={"shrink": 0.5})
plt.yticks(fontsize = 6)
plt.xticks(fontsize = 5)
plt.savefig('/home/wanglab/Desktop/inj.svg', bbox_inches = "tight")

#%%