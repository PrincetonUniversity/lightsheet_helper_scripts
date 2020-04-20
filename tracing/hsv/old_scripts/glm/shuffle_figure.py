#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 17:34:04 2019

@author: wanglab
"""

import seaborn as sns, numpy as np, os, matplotlib.pyplot as plt, pandas as pd
import pickle as pckl

""" 
dark blue is real data, red is what the model fit, the gray bars are the mean of all the 1000 shuffles
where the injection site was shuffled.
written in python 3.
bars on boxplots are standard deviation of the real data.
it is possible to do a t-test or something to see if we get significance that way, my guess is it is not necessary 
because we have the Wald t-test from the GLM also?
example of how to display shuffle using bar plots: https://www.nature.com/articles/nature23020/figures/4.
this paper does it side by side (shuffle next to data), but not sure if it is needed here again because we are not 
testing significance and just for visualization. it also helps with overcrowding to just overlay.
"""

#THALAMUS
#import data
pth = "/jukebox/wang/zahra/modeling/h129/thalamus/shuffle_figure_data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

#set dst 
dst = "/home/wanglab/Desktop"
   
#set alpha for gray transparency for shuffle
alpha = 0.3

#set colors for what should be what
actual = "darkblue"
fit = "firebrick"

#make boxplots to make fit and shuffle fit vs. actual
#change these if you need to ofc
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
#log was an idea but looks bad
df["shuffle_log"] = np.log10(data["fit_shuf"].mean(axis = 0).ravel())
df["fit_log"] = np.log10(df["count"].values)
df["brain"] = np.array(data["brains"]*18)
df["region"] = np.repeat(np.asarray(short_nuclei), 23)
df["inj"] = np.array([data["ak_pool"][idx] for idx in data["primary_pool"]]*18)

sns.set_style("white")
#removing simplex and lob iii, iv/v because n < 3
g = sns.FacetGrid(df[(df.inj != "Simplex") & (df.inj != "Lob. III, IV-V")], col="inj", height=5, aspect=1)

g.map(sns.barplot, "count", "region", facecolor=(1, 1, 1, 0), ci = "sd", linewidth=1, capsize = 0.2,
      edgecolor = actual) 
g.map(sns.swarmplot, "count", "region", color = actual, size = 3) #size determines radius of dot corresponding to data point
#g.map(sns.barplot, "fit", "region", facecolor=(1, 1, 1, 0), ci = None, capsize = 0.2, edgecolor = "firebrick") 
g.map(sns.swarmplot, "fit", "region", color = fit, size = 2) 
g.map(sns.barplot, "shuffle", "region", alpha = alpha, color = "gray", ci = None) 

g.set_xlabels("% of thalamic counts")
sns.despine(offset=10)

#save out    
plt.savefig(os.path.join(dst, "thal_shuffle_boxplot.pdf"), bbox_inches = "tight")

#NEOCORTEX
#import data
pth = "/jukebox/wang/zahra/modeling/h129/neocortex/shuffle_figure_data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

df = pd.DataFrame()
df["count"] = data["cell_counts_per_brain_p"].T.ravel()
df["fit"] = data["fit"].ravel()
df["shuffle"] = data["fit_shuf"].mean(axis = 0).ravel()
df["shuffle_log"] = np.log10(data["fit_shuf"].mean(axis = 0).ravel())
df["fit_log"] = np.log10(df["count"].values)
df["brain"] = np.array(data["brains"]*10)
df["region"] = np.repeat(np.asarray(data["regions"]), 33)
df["inj"] = np.array([data["ak_pool"][idx] for idx in data["primary_pool"]]*10)

#remove simplex, n < 3
df = df[df.inj != "Simplex"]

g = sns.FacetGrid(df[(df.region != "Somatomotor, Somatosensory")], 
    col="inj", height=3, aspect=1, sharex = True, sharey = True)

g.map(sns.barplot, "count", "region", facecolor=(1, 1, 1, 0), ci = "sd", linewidth=1, capsize = 0.2,
      edgecolor = actual) 
g.map(sns.swarmplot, "count", "region", color = actual, size = 3) 
#g.map(sns.barplot, "fit", "region", facecolor=(1, 1, 1, 0), ci = None, capsize = 0.2, edgecolor = "firebrick") 
g.map(sns.swarmplot, "fit", "region", color = fit, size = 2) 
g.map(sns.barplot, "shuffle", "region", alpha = alpha, color = "gray", ci = None) 

g.set_xlabels("% of neocortical counts")
sns.despine(offset=10)

#save out
plt.savefig(os.path.join(dst, "nc_shuffle_boxplot_except_sm.pdf"), bbox_inches = "tight")

#might need to plot sm separately as the x axis is too different - can merge these together in illustrator/ppt
g = sns.FacetGrid(df[(df.region == "Somatomotor, Somatosensory")], 
    col="inj", height=1.3, aspect=1.5, sharex = True, sharey = True)

g.map(sns.barplot, "count", "region", facecolor=(1, 1, 1, 0), ci = "sd", linewidth=1, capsize = 0.2,
      edgecolor = actual) 
g.map(sns.swarmplot, "count", "region", color = actual, size = 3)
#g.map(sns.barplot, "fit", "region", facecolor=(1, 1, 1, 0), ci = None, capsize = 0.2, edgecolor = "firebrick") 
g.map(sns.swarmplot, "fit", "region", color = fit, size = 2) 
g.map(sns.barplot, "shuffle", "region", alpha = alpha, color = "gray", ci = None) 
sns.despine(offset=10)

#save out
plt.savefig(os.path.join(dst, "nc_shuffle_boxplot_w_sm.pdf"), bbox_inches = "tight")