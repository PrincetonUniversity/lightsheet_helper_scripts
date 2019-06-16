#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 14:16:42 2019

@author: wanglab
"""

import seaborn as sns, numpy as np, os, matplotlib.pyplot as plt, pandas as pd
import pickle as pckl

#import data
pth = "/jukebox/wang/zahra/modeling/h129/thalamus/shuffle_figure_data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

#set dst 
dst = "/home/wanglab/Desktop"
    
#make boxplots to make fit and shuffle fit vs. actual?
short_nuclei = ['Parafascicular n.',
 'P. complex',
 'P. triangular n.',
 'LP',
 'L. habenula',
 'LD',
 'CL',
 'PV',
 'Reuniens',
 'MD',
 'V. of lat. gen. complex',
 'VPL',
 'VPM',
 'Submedial n.',
 'RTN',
 'VM',
 'AV',
 'VA-L']

df = pd.DataFrame()
df["count"] = data["cell_counts_per_brain_p"].T.ravel()
df["fit"] = data["fit"].ravel()
df["shuffle"] = data["fit_shuf"].mean(axis = 0).ravel()
df["shuffle_log"] = np.log10(data["fit_shuf"].mean(axis = 0).ravel())
df["fit_log"] = np.log10(df["count"].values)
df["brain"] = np.array(data["brains"]*18)
df["region"] = np.repeat(np.asarray(short_nuclei), 23)
df["inj"] = np.array([data["ak_pool"][idx] for idx in data["primary_pool"]]*18)

sns.set_style("white")

g = sns.FacetGrid(df, col="inj", height=5, aspect=.5)

g.map(sns.barplot, "fit", "region", facecolor=(1, 1, 1, 0), ci = "sd", capsize = 0.2, edgecolor = "darkblue") 
g.map(sns.swarmplot, "fit", "region", color = "darkblue", size = 3) 
g.map(sns.barplot, "shuffle", "region", alpha = 0.4, color = "gray", ci = None) 

g.set_xlabels("% of thalamic counts")

sns.despine(offset=10)
    
plt.savefig(os.path.join(dst, "boxplot.pdf"), bbox_inches = "tight")