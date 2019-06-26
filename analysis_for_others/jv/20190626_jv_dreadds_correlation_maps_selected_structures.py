#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 18:40:33 2019

@author: wanglab
"""

import pandas as pd, os, matplotlib as mpl, matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

"""
Calculates a Spearman rank-order correlation coefficient and the p-value to test for non-correlation.

The Spearman correlation is a nonparametric measure of the monotonicity of the relationship between two datasets. 
Unlike the Pearson correlation, the Spearman correlation does not assume that both datasets are normally distributed. 
Like other correlation coefficients, this one varies between -1 and +1 with 0 implying no correlation. 
Correlations of -1 or +1 imply an exact monotonic relationship. Positive correlations imply that as x increases, so does y. 
Negative correlations imply that as x increases, y decreases.

The p-value roughly indicates the probability of an uncorrelated system producing datasets that have a Spearman correlation 
at least as extreme as the one computed from these datasets. The p-values are not entirely reliable but are probably reasonable 
for datasets larger than 500 or so.
"""

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

### manipulating dataframe
#################################################################CHANGE PATHS##################################################################
#set destination for figs
dst = "/home/wanglab/Desktop"

#set source
src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis"
#import counts
pth = os.path.join(src, "select_structures_percent_counts.csv")

cts = pd.read_csv(pth)
#wtf - rename back to original
cts["animal"] = cts["Unnamed: 0"]
cts = cts.drop(columns = ["Unnamed: 0"])

#import behavior & inj metrics
bh_pth = os.path.join(src, "ymaze.csv")
bh = pd.read_csv(bh_pth)

#sort by animal in counts csv
cts = cts.sort_values(by = ["animal"])

#remove homecage control for now from counts
cts = cts[cts.group != "homecage_control"]

#join both dataframes for count and behavior?
df = pd.concat([cts, bh], axis = 1)

#drop unneeded
df = df.drop(columns = ["Mouse ID", "Condition", "rev", "Region"])

#only DREADDS
df = df[df.group == "DREADDs"]

#corr
#get structures
structs = [xx[0] for xx in pd.read_csv(os.path.join(src, "structures.csv")).values]
a = df.lobvi.values #this is the things you are correlating each structure's % counts to, 
#theoretically you can change this to a behavior metric to see how it changes the coefficients....

#init empty df
dfcorr = pd.DataFrame()
dfcorr["name"] = structs
for struct in structs:
    b = df[struct].values
    rho = spearmanr(a, b)[0] #correlation
    pval = spearmanr(a, b)[1] #pvalue
    dfcorr.loc[(dfcorr.name == struct), "corr_coeff"] = rho
    dfcorr.loc[(dfcorr.name == struct), "corr_pval"] = pval

dfcorr.to_csv(os.path.join(dst, "spearman_coeff_select_structures_dreadds.csv"))

#%%
#plotting - NOT complete
sig_structs = [xx for xx in dfcorr.name.values if dfcorr.loc[(dfcorr.name == xx), "corr_pval"].values[0] < 0.05]

counts = np.asarray([df[xx].values for xx in sig_structs])

regions = np.asarray(sig_structs)

coeff = np.asarray([dfcorr.loc[(dfcorr.name == xx), "corr_coeff"].values[0] for xx in dfcorr.name.values if dfcorr.loc[(dfcorr.name == xx), "corr_pval"].values[0] < 0.05])

brains = df.animal.values

#formatting
fig = plt.figure(figsize=(3,6))
ax = fig.add_axes([.4,.1,.5,.7])

show = counts

vmin = 0
vmax = 3#.6
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,3,7)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("% of total counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

#value annotations
for ri,row in enumerate(coeff):
    for ci,col in enumerate(row):
        pass
        if counts[ < 25:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")
        
        
cb.ax.set_visible(True)
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}".format(br) for br in brains], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="xx-small")
plt.savefig(os.path.join(dst,"cfos_percent_counts_1.pdf"), bbox_inches = "tight")