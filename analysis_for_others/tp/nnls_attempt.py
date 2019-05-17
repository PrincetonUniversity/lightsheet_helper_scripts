#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:26:18 2019

@author: wanglab
"""

from scipy.optimize import nnls

mat = []
mat_shuf = []
rnrms = []
rnrms_shuf = []
for itera in range(1000):
    print(itera)
    if itera == 0:
        shuffle = False
        inj = X.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        rnrms_shuf.append([])
        inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:]
    for count, region in zip(Y.T, regions):
        inj_ = inj[~np.isnan(count)]
        count = count[~np.isnan(count)]

        # intercept:
        inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
        
        w, resnorm = nnls(inj_, count)

        if not shuffle:
            mat.append(w[:-1])
            rnrms.append(resnorm)
        elif shuffle:
            mat_shuf[-1].append(w[:-1])
            rnrms_shuf[-1].append(resnorm)
        
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
rnrms_shuf = np.array(rnrms_shuf)
#%%
## display
fig = plt.figure()#figsize=(4,7))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat # NOTE abs
amax = 4.7 #computed from np.max(np.abs(mat))
#amax = np.max(np.abs(mat))
#vmin = -amax
vmin = np.min(mat)-0.5
vmax = np.max(mat)+0.5
cmap = plt.cm.RdBu_r
cmap.set_under('w')
#colormap
#cmap = plt.cm.bwr
# discrete colorbar details
bounds = np.linspace(-5,5,11)
#bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", shrink=0.5, aspect=10)
cb.set_label("Weight", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
# xticks
ax.set_xticks(np.arange(len(ak))+.5)

#remaking labeles so it doesn't look squished
lbls = np.asarray(ak)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
#The adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases
# only if the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model 
# by less than expected by chance. The adjusted R-squared can be negative, but it’s usually not.  It is always lower than the R-squared.
#ax.set_yticklabels(["{}\nr2adj={:0.2f}".format(bi,ar) for ar,bi in zip(ars,regions)], fontsize="xx-small")
ax.set_yticklabels(["{}".format(bi) for bi in regions], fontsize="xx-small")
plt.savefig("/home/wanglab/Desktop/weights.svg", bbox_inches = "tight")
