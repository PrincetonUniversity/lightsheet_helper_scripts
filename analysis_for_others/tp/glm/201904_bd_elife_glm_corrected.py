#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 15:09:22 2019

@author: wanglab
"""

## load data

from scipy.io import loadmat
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
from scipy.stats import spearmanr

datafile = '/home/wanglab/Desktop/data.mat'
data = loadmat(datafile)

## preprocess data

age = data['sample_mice_age']
age = np.array([str(i[0].squeeze()) for i in age])
age[age=='Juven'] = 0
age[age=='Adult'] = 1
age = age.astype(int)

behav = data['behavior']
bnames = data['behavior_metric_names']
bnames = np.array([str(i.squeeze()) for i in bnames[0]])

expr = data['sample_signal']
expr_raw = np.array([e[0] for e in expr])
expr_raw = expr_raw.astype(bool)

atlas = data['atlas_rois'] # y, x, d
atlas = np.array([a[0] for a in atlas]) # lobule, xyz
atlas = atlas.astype(bool)
atlas_keys = data['atlas_roi_names']
atlas_keys = np.array([str(i.squeeze()) for i in atlas_keys[0]])

expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for lob in atlas] for mouse in expr_raw])
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for lob in atlas] for mouse in expr_raw])
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

# categorizing (not necessary for most things so far)
of_interest = np.array([1,3,5,6])
primary_among_4_of_interest_coded = np.argmax(expr_all_as_frac_of_lob[:,of_interest], axis=1)
primary_among_4_of_interest = of_interest[primary_among_4_of_interest_coded]

## compute stat for paper text
for lid,lstr in zip([5,6,1,3],['lob6','lob7','crus1','crus2']):
    grp = expr_all_as_frac_of_lob[primary==lid][:,lid]
    mean = grp.mean()
    std = grp.std(ddof=1)
    print('{}\t{:0.2f} ± {:0.2f}'.format(lstr, mean, std))

frac2 = np.array([ms[scd] for ms,scd in zip(expr_all_as_frac_of_lob, secondary)])
print('\nFraction of lobule for secondary lobules: {:0.2f} ± {:0.2f}'.format(frac2.mean(), frac2.std(ddof=1)))

print('\nGropuing by primary lobule, mean fraction of secondary lobule covered:')
for p in np.unique(primary):
    scd = secondary[primary==p]
    mc = expr_all_as_frac_of_lob[primary==p]
    fracs = [mi[si] for mi,si in zip(mc,scd)]
    print('primary =',atlas_keys[p])
    print('\tfraction of secondary = {:0.2f}'.format(np.mean(fracs)))

## split lobules
lob_i = [1,3,5,6,9]
lob_names = ['R_CrusI', 'R_CrusII', 'LobVI', 'LobVII', 'R_PM']
ml_splits = [1, 1, 1, 1, 1]

include = [0,1,2,3] # index of 3 lists directly above this
lob_i = np.array(lob_i)[include]
lob_names = np.array(lob_names)[include]
ml_splits = np.array(ml_splits)[include]

lobs = atlas[lob_i]
ak = []
atlas_mod = []
for name,split,lob in zip(lob_names, ml_splits, lobs):
    lobpos = np.where(lob)
    ml_pos = lobpos[1]
    ml_lims = np.round(np.linspace(min(ml_pos), max(ml_pos)+1, split+1)).astype(int)
    for spl in range(split):
        i0,i1 = ml_lims[spl:spl+2]
        coords = tuple([l[(ml_pos>=i0) & (ml_pos<i1)] for l in lobpos])
        sublob = np.zeros_like(lob)
        sublob[coords] = 1

        atlas_mod.append(sublob)
        if split > 1:
            ak.append(name + '_' + str(spl))
        else:
            ak.append(name)
atlas_mod = np.array(atlas_mod)

## compute expressions

expr = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for lob in atlas_mod] for mouse in expr_raw])

## filters
include_age = 1# 0=juv, 1=adult
behav_exclude = [ 
                  'YM_AcqAbility',
                  #'YM_AcqInitialLR',
                  #'YM_AcqSecondaryLR',
                  #'YM_EarlyRevAbility',
                  #'YM_EarlyRevInitialLR',
                  #'YM_EarlyRevSecondaryLR',
                  'YM_LateRevAbility',
                  'YM_LateRevInitialLR',
                  'YM_Velocity',
                  'EPM_Velocity',
                  ]

e = expr.copy()
ex__ = e[age==include_age] # mice x voxel

b = behav.copy()
bh = b[age==include_age] # mice x behavior

prim = primary.copy()
prim = prim[age==include_age]

bhn = bnames.copy()
binclude = np.array([bi not in behav_exclude for bi in bhn])
bh = bh[:,binclude]
bhn = bhn[binclude]

#np.random.shuffle(bh)

##  glm
mat = []
pmat = []
mat_shuf = []
p_shuf = []
ars = []
for itera in range(100):
    if itera == 0:
        shuffle = False
        ex = ex__.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        ex = ex__[np.random.choice(np.arange(len(ex)), replace=False, size=len(ex)),:]
    for b,bname in zip(bh.T, bhn):
        ex_ = ex[~np.isnan(b)]
        b = b[~np.isnan(b)]

        # intercept:
        ex_ = np.concatenate([ex_, np.ones(ex_.shape[0])[:,None]], axis=1)

        glm = sm.OLS(b, ex_)
        #glm_fam = sm.families.Gaussian(sm.genmod.families.links.identity)
        #glm = sm.GLM(b, ex_, family=glm_fam)
        res = glm.fit()

        coef = res.params[:-1]
        se = res.bse[:-1]
        pvals = res.pvalues[:-1]

        val = coef/se

        if not shuffle:
            mat.append(val)
            pmat.append(pvals)
            ars.append(res.rsquared_adj)
        elif shuffle:
            mat_shuf[-1].append(val)
            p_shuf[-1].append(pvals)
        
        # inspect residuals
        """
        if not shuffle:
            pl.clf()
            pl.scatter(res.fittedvalues, res.resid)
            pl.hlines(0, res.fittedvalues.min(), res.fittedvalues.max())
            pl.title(bname)
            pl.xlabel('Fit-vals')
            pl.ylabel('Residuals')
            pl.savefig('figs/resid_inspection-{}.png'.format(bname))
        """

mat = np.array(mat) # behavior x voxel
mat_shuf = np.array(mat_shuf) # behavior x voxel
pmat = np.array(pmat) # behavior x voxel
p_shuf = np.array(p_shuf)
ars = np.array(ars)

#ordr = np.argsort(np.max(np.abs(mat),axis=1))
ordr = np.arange(len(mat))
mat = mat[ordr,:]
mat_shuf = mat_shuf[:,ordr,:]
p_shuf = p_shuf[:,ordr,:]
pmat = pmat[ordr,:]
bn = bhn[ordr]
#%%
## display
fig = pl.figure(figsize=(4,7))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = mat # NOTE abs
amax = 4.47 #computed from juvenile, so colors match for both #np.max(np.abs(mat))
#amax = np.max(np.abs(mat))
#vmin = -amax
vmin = -3
vmax = amax
#cmap = pl.cm.RdBu_r
cmap = pl.cm.Reds

# discrete colorbar details
#bounds = np.linspace(-5,5,11)
bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label='Weight / SE', shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i', shrink=0.5, aspect=10)
cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%0.1f', shrink=0.5, aspect=10)
cb.set_label('| Weight / SE |', fontsize='x-small', labelpad=3)
cb.ax.tick_params(labelsize='x-small')

if include_age == 1:
    cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        ax.text(ci+.5, ri+.5, '{:0.1f}'.format(col), color='k', fontsize=8, ha='center', va='center')

# signif
sig = pmat < .05#/np.size(pmat)
null = (p_shuf < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
#    ax.text(x+.5, y+.5, 'x', fontsize=15, ha='center', va='center', transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) X's are expected by chance if no real effect exists".format(nullmean, nullstd), ha='center', va='center', fontsize='small', transform=ax.transAxes)

# aesthetics
# xticks
ax.set_xticks(np.arange(len(ak))+.5)
ax.set_xticklabels(ak, rotation=30, fontsize='x-small', ha='center')
# yticks
ax.set_yticks(np.arange(len(bhn))+.5)
#ax.set_yticklabels(['{}\nr2adj={:0.2f}'.format(bi,ar) for ar,bi in zip(ars,bn)], fontsize='x-small')
ax.set_yticklabels(['{}'.format(bi) for ar,bi in zip(ars,bn)], fontsize='x-small')

#%%
## specific relevant scatter plots

# juvenile
if include_age == 0:
    #bchoice = 'SC_NoveltySeeking'
    #bchoice = 'SC_BSDistance'
    bchoices = bn.copy()

    echoices = ['LobVI','R_CrusI','R_CrusII','LobVII']
    primids = [5,1,3,6]

    choice_idx = 1 # default 1
    echoice = echoices[choice_idx]
    col = (prim==primids[choice_idx]).astype(int) # should match echoice

elif include_age == 1:
    #bchoice = 'YM_EarlyRevInitialLR'
    #bchoice = 'GR_GroomingRatio'
    #bchoice = 'YM_Distance'
    #bchoice = 'EPM_AntiHesitation'
    #bchoice = 'YM_AcqInitialLR'
    bchoices = bn.copy()

    echoices = ['LobVI','R_CrusI','R_CrusII','LobVII']
    primids = [5,1,3,6]
    
    choice_idx = 0 # from list directly above this, default 0
    echoice = echoices[choice_idx]
    col = (prim==primids[choice_idx]).astype(int) # should match echoice


fig,ax = pl.subplots(1,1,figsize=(5,5), gridspec_kw=dict(bottom=.2, left=.2, right=.9, top=.9))
for bchoice in bchoices:
    ax.clear()

    bi = int(np.where(bn == bchoice)[0].squeeze())
    ei = int(np.where(np.array(ak) == echoice)[0].squeeze())
    ex_vec = ex__[:,ei]
    b_vec = bh[:,bi]

    mfc = ['none','k']

    for mfci,idx in zip(mfc, [0,1]):
        ax.scatter(ex_vec[col==idx], b_vec[col==idx], color='k', facecolor=mfci, linewidth=.2, s=160)
        ax.set_xlabel('Fraction of {} labelled'.format(echoice))
        ax.set_ylabel('{}'.format(bchoice))
    pl.savefig('/home/wanglab/Desktop/scatter_{}_{}_{}.pdf'.format({0:'juv',1:'adult'}[include_age], bchoice,echoice))

## save out all inspection scatterplots
pl.close('all')

for lob,lobname in zip(ex.T, atlas_keys):
    print(lobname)
    fig,axs = pl.subplots(3,7)
    for bname,beh,ax in zip(bhn, bh.T, axs.ravel()):
        ax.scatter(lob,beh)
        ax.set_title(bname, fontsize='xx-small')
        ax.set_xticks([])
        ax.set_yticks([])
    pl.figtext(.5, .05, 'DREADD expression', fontsize='small')
    pl.figtext(.1, .5, 'Behavior metric', fontsize='small', rotation=90)
    fig.savefig('/home/wanglab/Desktop/bylob-'+lobname.replace('/','')+'.png')
    pl.close(fig)

for bname,beh in zip(bhn, bh.T):
    print(bname)
    fig,axs = pl.subplots(3,4)
    for lob,lobname,ax in zip(ex.T, atlas_keys, axs.ravel()):
        ax.scatter(lob,beh)
        ax.set_title(lobname, fontsize='xx-small')
        ax.set_xticks([])
        ax.set_yticks([])
    pl.figtext(.5, .05, 'DREADD expression', fontsize='small')
    pl.figtext(.1, .5, 'Behavior metric', fontsize='small', rotation=90)
    fig.savefig('/home/wanglab/Desktop/bybehav-'+bname+'.png')
    pl.close(fig)
## correlation matrices
lname = 'Crus 1 (Right)'
lname = 'Lobule VI'
rho0 = spearmanr(behav[(atlas_keys[primary]==lname) & (age==0)])[0]
rho1 = spearmanr(behav[(atlas_keys[primary]==lname) & (age==1)])[0]
rho0[np.isnan(rho0)] = 0
rho1[np.isnan(rho1)] = 0
fig,ax = pl.subplots(1,1)
ax.pcolor(np.abs(rho0-rho1))
##