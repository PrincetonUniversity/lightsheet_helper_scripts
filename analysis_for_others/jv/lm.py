#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 15:28:35 2019

@author: wanglab
"""

import pandas as pd, numpy as np, os
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns 
from scipy.stats import spearmanr

csv_pth = r"C:\Users\Zahra\Desktop\select_structures_percent_counts_for_visualization.csv"#r"C:\Users\Zahra\Desktop\select_structures_percent_counts_for_visualization.csv"

cfos = pd.read_csv(csv_pth)
cfos["animal"] = [int(xx[2:]) for xx in cfos.animal.values]
cfos = cfos.sort_values(by = ["animal"])
#look at certain groups
cfos = cfos[cfos.group != "homecage_control"]

#look at nc first?
#ncs = ["Anterior cingulate area", "Prelimbic area", "Infralimbic area",
#       "Orbital area"]

#striatum
#ncs = ["Striatum dorsal region", "Striatum ventral region","Striatum-like amygdalar nuclei", "Pallidum"]

#amygdala
#ncs = ["Lateral amygdalar nucleus", "Basolateral amygdalar nucleus",
#       "Basomedial amygdalar nucleus", "Posterior amygdalar nucleus"]
##
#thalamus & hypothalamus
#ncs = ["Ventral group of the dorsal thalamus",
#       "Subparafascicular nucleus",
#       "Subparafascicular area", "Peripeduncular nucleus",
#       "Geniculate group, dorsal thalamus",
#ncs = ["Lateral group of the dorsal thalamus",
#       "Anterior group of the dorsal thalamus",
#       "Medial group of the dorsal thalamus",
#       "Midline group of the dorsal thalamus",
#       "Intralaminar nuclei of the dorsal thalamus",
#       "Reticular nucleus of the thalamus"]
#       "Geniculate group, ventral thalamus"]

##midbrain
#ncs = ["Midbrain, sensory related", "Midbrain, motor related",
#       "Midbrain, behavioral state related"]

##cortical subplate
#ncs = ["Cortical amygdalar area","Hippocampal region", "Retrohippocampal region"]
#
ncs = [xx for xx in list(cfos.columns) if not xx == "animal" and not xx == "group" and cfos[xx][8] > 0]

#x
counts = np.asarray([cfos[nc].values for nc in ncs]); counts = counts.T

##combine counts for orbitofrontal cortex
#counts_pl = np.asarray([[brain[0], brain[1] + brain[2] + brain[3] + brain[4], brain[5], brain[6], brain[7]] for brain in counts])
# 
##rename regions
#ncs_pl = ["Somatosensory areas", 
#       "Infralimbic, Prelimbic, Anterior cingulate, Orbital",
#       "Agranular insular area", "Retrosplenial area",
#       "Posterior parietal association areas"]

ymaze_pth = r"C:\Users\Zahra\Desktop\ymaze.csv"#r"C:\Users\Zahra\Desktop\ymaze.csv"

ym = pd.read_csv(ymaze_pth)

#look at certain groups
#ym = ym[ym.Condition != "no rev"]

mts = ["Initial Reversal",
       "Multisession Reversal ", "Final Reversal", "Forced Reversal"]
#y
ym_mt = np.asarray([ym[mt].values for mt in mts])

#variables
X = counts
Y = ym_mt

#%%

#set destination
dst = r"C:\Users\Zahra\Desktop"

X = counts
X = np.rot90(X, 1)

ncs = ncs[::-1]
cmap = plt.cm.viridis

vmin = 0
vmax = 8
cmap.set_under('w')
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize = (6, 7))
pc = ax.pcolor(X, cmap=cmap, vmin=vmin, vmax=vmax)
ax.set_yticks(np.arange(len(ncs)))
ax.set_yticklabels(["{}".format(ak) for ak in ncs], fontsize=7, va = "bottom", ha = "right")
ax.set_xticks(np.arange(len(cfos["animal"].values)))
ax.set_xticklabels(["{}".format(brain) for brain in cfos["animal"]], rotation=30, fontsize=5, ha="right")
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.3, aspect=10)
cb.ax.tick_params(labelsize="x-small")
cb.set_label("% of total cell counts", fontsize="x-small", labelpad=3)

import matplotlib.patches as mpatches

green_patch = mpatches.Patch(color="navy", label="DREADDs")
gold_patch = mpatches.Patch(color="darkgreen", label="CNO reversal")
blue_patch = mpatches.Patch(color="darkred", label="CNO no reversal")

legend = plt.legend(handles=[green_patch, gold_patch, blue_patch], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

frame = legend.get_frame()
frame.set_facecolor('white')

plt.savefig(os.path.join(dst, "counts.svg"), bbox_inches = "tight")

#%%
plt.figure(figsize = (6, 0.6))
sns.set(font_scale=0.5)
mts = ["Initial Reversal","Final Reversal", "Forced Reversal"]
#y
ym_mt = np.asarray([ym[mt].values for mt in mts])

#accompanying heatmap
mp = pd.DataFrame(ym_mt.T)
mp = mp.set_index(cfos.animal.values)
mp.columns = mts

sns.heatmap(mp.T, cmap = "Reds_r")
plt.savefig(os.path.join(dst, "behav.svg"), bbox_inches = "tight")

#%%
plt.figure(figsize = (6, 0.3))

#plot groups
cfos.loc[(cfos.group == "CNO_control_no_reversal"), "groups"] = 1
cfos.loc[(cfos.group == "CNO_control_reversal"), "groups"] = 2
cfos.loc[(cfos.group == "DREADDs"), "groups"] = 3

sns.heatmap(pd.DataFrame(cfos.groups.values).T, cmap = sns.hls_palette(3, l=.3, s=.8))

plt.savefig(os.path.join(dst, "groups.svg"), bbox_inches = "tight")

#%%
#from badura et al.
mat = []
pmat = []
mat_shuf = []
p_shuf = []
ars = []
rs = []
for itera in range(1000):
    if itera%100 == 0: print(itera)
    if itera == 0:
        shuffle = False
        count = X.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        count = X[np.random.choice(np.arange(len(count)), replace=False, size=len(count)),:]
    for ymaze, metric in zip(Y, mts):

        # intercept:
        count_ = np.concatenate([count, np.ones(count.shape[0])[:,None]*1], axis=1)
        
        glm = sm.OLS(ymaze, count_)
#        glm = sm.GLM(count, inj_, family=sm.families.Poisson())
        res = glm.fit()
        
        coef = res.params[:-1]
        se = res.bse[:-1]
        pvals = res.pvalues[:-1] 

        val = coef/se

        if not shuffle:
            mat.append(val)
            pmat.append(pvals)
            ars.append(res.rsquared_adj)
            rs.append(res.rsquared)
        elif shuffle:
            mat_shuf[-1].append(val)
            p_shuf[-1].append(pvals)
        
        # inspect residuals
        """
        if not shuffle:
            plt.clf()
            plt.scatter(res.fittedvalues, res.resid)
            plt.hlines(0, res.fittedvalues.min(), res.fittedvalues.max())
            plt.title(region)
            plt.xlabel("Fit-vals")
            plt.ylabel("Residuals")
            plt.savefig(os.path.join(dst, "resid_inspection-{}.png").format(region))
        """
        

mat = np.array(mat) # region x inj
mat_shuf = np.array(mat_shuf) # region x inj
pmat = np.array(pmat) # region x inj
p_shuf = np.array(p_shuf)
ars = np.array(ars)
rs = np.array(rs)

#%%
## display
fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat # NOTE abs

vmin = -3
vmax = 3
cmap = plt.cm.RdBu_r
cmap.set_under("steelblue")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,vmax-vmin+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", shrink=0.5, aspect=10)
cb.set_label("Weight / SE", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# signif
sig = pmat < .05#/np.size(pmat)
#p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05", ha="center", va="center", fontsize="small", transform=ax.transAxes)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
# xticksjt -t monokai -m 200
ax.set_xticks(np.arange(len(ncs))+.5)

#remaking labeles so it doesn"t look squished
ax.set_xticklabels(["{}".format(nc) for nc in ncs], rotation=20, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(mts))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(mt) for mt in mts], fontsize="xx-small")
plt.savefig("/home/wanglab/Desktop/weights.pdf", bbox_inches = "tight")
