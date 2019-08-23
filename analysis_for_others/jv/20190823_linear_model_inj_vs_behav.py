#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 16:21:15 2019

@author: wanglab
"""

import pandas as pd, numpy as np, statsmodels.api as sm, matplotlib.pyplot as plt, matplotlib as mpl, os

df_pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/more_selected_structures/ymaze.csv"

df = pd.read_csv(df_pth)

df = df[df.Region == "Lobule VI"]
inj = np.asarray([df.lobvi.values, df.lobvii.values, df.crus1.values, df.crus2.values]).T
bh = np.asarray([df["InitialReversal"], df["Multisession Reversal"], df["FinalReversal"], df["ForcedReversal"]]).T
X = inj
Y = bh
regions = ["Lob. VI", "Lob. VII", "Crus I", "Crus II"]
bhn = ["Initial Reversal", "Multisession Reversal", "Final Reversal", "Forced Reversal"]

dst = "/home/wanglab/Desktop"
c_mat = []
mat = []
pmat = []
mat_shuf = []
p_shuf = []
fit = []
fit_shuf = []

for itera in range(1000):
    if itera%100 == 0: print(itera)
    if itera == 0:
        shuffle = False
        inj = X.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        fit_shuf.append([])
        inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:]
    for count, region in zip(Y.T, regions):
        inj_ = inj[~np.isnan(count)]
        count = count[~np.isnan(count)]

        # intercept:
        inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
        
        glm = sm.OLS(count, inj_)
#        glm = sm.GLM(count, inj_, family=sm.families.Poisson())
        res = glm.fit()
        
        coef = res.params[:-1]
        se = res.bse[:-1]
        pvals = res.pvalues[:-1] 

        val = coef/se

        if not shuffle:
            c_mat.append(coef)
            mat.append(val)
            pmat.append(pvals)
            fit.append(res.fittedvalues)
            
        elif shuffle:
            mat_shuf[-1].append(val)
            p_shuf[-1].append(pvals)
            fit_shuf[-1].append(res.fittedvalues)
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
        
c_mat = np.array(c_mat)
mat = np.array(mat) # region x inj
mat_shuf = np.array(mat_shuf) # region x inj
pmat = np.array(pmat) # region x inj
p_shuf = np.array(p_shuf)
fit = np.array(fit)
fit_shuf = np.array(fit_shuf)

#%%
## display
fig = plt.figure(figsize=(5,2.5))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat 

vmin = -1
vmax = 2
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,9)
#bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds,
                  boundaries=bounds, format="%0.1f", shrink=0.5, aspect=10)
cb.set_label("Weight / SE", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(mat):
    for ci,col in enumerate(row):
        if col > 0:
            ax.text(ci+.5, ri+.5, "{:0.2f}".format(col), color="w", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.2f}".format(col), color="k", ha="center", va="center", fontsize="x-small")
        
# signif
sig = pmat < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)

ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, 
        nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
# xticksjt -t monokai -m 200
ax.set_xticks(np.arange(len(regions))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(regions)
ax.set_xticklabels(lbls, rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(bhn))+.5)
ax.set_yticklabels(bhn, fontsize="xx-small")
plt.savefig(os.path.join(dst, "ymaze_lm.pdf"), bbox_inches = "tight")

