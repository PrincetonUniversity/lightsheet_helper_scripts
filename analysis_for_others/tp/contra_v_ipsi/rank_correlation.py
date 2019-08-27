#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 17:21:19 2019

@author: wanglab
"""

import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, matplotlib as mpl, statsmodels.api as sm, matplotlib.ticker as ticker

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#import data
main_data_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_contra_ipsi_counts_densities.p"
data = pckl.load(open(main_data_pth, "rb"), encoding = "latin1")

#injection site analysis
inj_pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
inj_dct = pckl.load(open(inj_pth, "rb"), encoding = "latin1")

#inj volumes
inj_vol_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_inj_vol.p" 
inj_vol_dct = pckl.load(open(inj_vol_pth, "rb"), encoding = "latin1")
inj_vol = inj_vol_dct["inj_vol"]
iv = []
for k,v in inj_vol.items():
    iv.append(v)

#set save destination
sv_dst = "/home/wanglab/Desktop"

cell_counts_per_brain_left = data["cell_counts_per_brain_left"]
cell_counts_per_brain_right = data["cell_counts_per_brain_right"]
density_per_brain_left = data["density_per_brain_left"]
density_per_brain_right = data["density_per_brain_right"]
lr_dist = data["lr_dist"]
nc_areas = data["nc_areas"] #gives order of nc areas also

brains = inj_dct["brainnames"]
primary_pool = inj_dct["primary_pool"]
ak_pool = inj_dct["cb_regions_pool"]
inj = inj_dct["expr_all_as_frac_of_inj_pool"]
vols = [xx/1e9 for xx in iv]

#for rank correlation calculation
total_cell_counts_per_brain = np.sum(cell_counts_per_brain_left+cell_counts_per_brain_right, axis=1)

#preprocessing
nc_left_counts = np.sum(cell_counts_per_brain_left, axis=1)
nc_right_counts = np.sum(cell_counts_per_brain_right, axis=1)
nc_density_left = density_per_brain_left
nc_density_right = density_per_brain_right

lrv = list(lr_dist.values())
lr_brains = list(lr_dist.keys())

#dct is just for my sanity, so im not mixing up brains
_ccontra = []; _cipsi = []; _dcontra = []; _dipsi = []
for i in range(len(lr_brains)):
    if lrv[i] > 0: #right
        #counts
        _ccontra.append(nc_left_counts[i])
        _cipsi.append(nc_right_counts[i])
        #density
        _dcontra.append(nc_density_left[i])
        _dipsi.append(nc_density_right[i])
    elif lrv[i] < 0: #left
        #counts
        _ccontra.append(nc_right_counts[i])
        _cipsi.append(nc_left_counts[i])
        #density
        _dcontra.append(nc_density_right[i])
        _dipsi.append(nc_density_left[i])


_ccontra = np.asarray(_ccontra).T; _dcontra = np.asarray(_dcontra).T
_cipsi = np.asarray(_cipsi).T; _dipsi = np.asarray(_dipsi).T
_dratio = np.asarray([_dcontra[i]/_dipsi[i] for i in range(len(_dcontra))])
_cratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
#make into one
_dist = np.asarray(list(lr_dist.values()))

#injection site analysis
pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

brains = data["brainnames"]
primary_pool = data["primary_pool"]
ak_pool = data["cb_regions_pool"]
inj = data["expr_all_as_frac_of_inj_pool"]
 
_inj = np.asarray([inj[i] for i in range(len(inj)) if brains[i] in lr_brains])
_primary_pool = np.asarray([primary_pool[i] for i in range(len(primary_pool)) if brains[i] in lr_brains])
#%%
fig = plt.figure()
ax = fig.add_axes([.15,.15,.5,.8])

X = total_cell_counts_per_brain
Y = _cratio.T

ax.scatter(X, Y, color = "k", s = 10)

ax.set_xlabel("Total neocortical cell counts")
ax.set_ylabel("Contra/Ipsi count ratio")
ax.axhline(y=1.2, color = "grey")
plt.savefig(os.path.join(sv_dst, "nc_contra_ipsi_raw_scatter.pdf"), dpi=300)

fig = plt.figure(figsize=(12,7))
ax = fig.add_axes([.15,.15,.5,.8])

#calculate ranks
X = np.argsort(np.sort(total_cell_counts_per_brain))+1
Y = _cratio.T[np.argsort(np.sort(total_cell_counts_per_brain))].argsort()+1 #check if works

results = sm.OLS(Y,sm.add_constant(X)).fit()

m = results.params[1]
r2 = results.rsquared
b = results.params[0]

textstr = "\n".join((
    "slope: {:0.2f}".format(m),
    "$R^2$: {:0.2f}".format(r2)))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

ax.scatter(X, Y, color = "red", s = 10)
ax.plot(m*X + b, color = "k", linestyle = "--")

ytick_spacing = 1; xtick_spacing = 1
ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
ax.set_xlim([0,34])
ax.set_xlabel("Neocortical counts rank order")
ax.set_ylabel("Contra/Ipsi count ratio rank order")
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment="top", bbox=props)
plt.savefig(os.path.join(sv_dst, "nc_contra_ipsi_rank.pdf"), dpi=300)

#%%
#import data
main_data_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_contra_ipsi_counts_densities.p"
data = pckl.load(open(main_data_pth, "rb"), encoding = "latin1")

#injection site analysis
inj_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data.p"
inj_dct = pckl.load(open(inj_pth, "rb"), encoding = "latin1")

#inj volumes
inj_vol_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_inj_vol.p" 
inj_vol_dct = pckl.load(open(inj_vol_pth, "rb"), encoding = "latin1")
inj_vol = inj_vol_dct["inj_vol"]
iv = []
for k,v in inj_vol.items():
    iv.append(v)
vols = [xx/1e9 for xx in iv]

#set dst 
sv_dst = "/home/wanglab/Desktop"

#mask unwanted brains - dropping fast spreading ones and a low count one
curated_brains = [True, False, False, True, True, True, False, True, True, True, True, True,
                   True, True, True, True, True, True, True, True, True, True, True]

cell_counts_per_brain_left = data["cell_counts_per_brain_left"][curated_brains]
cell_counts_per_brain_right = data["cell_counts_per_brain_right"][curated_brains]
density_per_brain_left = data["density_per_brain_left"][curated_brains]
density_per_brain_right = data["density_per_brain_right"][curated_brains]
lr_dist = data["lr_dist"]

brains = np.array(inj_dct["brainnames"])[curated_brains]
primary_pool = inj_dct["primary_pool"][curated_brains]
ak_pool = inj_dct["cb_regions_pool"]
inj = inj_dct["expr_all_as_frac_of_inj_pool"][curated_brains]

#for rank correlation calculation
total_cell_counts_per_brain = np.sum(cell_counts_per_brain_left+cell_counts_per_brain_right, axis=1)

#-------------------------------------------------------------------------------------------------------------------------------------
#preprocessing
thal_left_counts = np.sum(cell_counts_per_brain_left, axis=1)
thal_right_counts = np.sum(cell_counts_per_brain_right, axis=1)
thal_density_left = density_per_brain_left
thal_density_right = density_per_brain_right

lrv = np.array(list(lr_dist.values()))[curated_brains]
lr_brains = np.array(list(lr_dist.keys()))[curated_brains]

_ccontra = []; _cipsi = []; _dcontra = []; _dipsi = []
for i in range(len(lr_brains)):
    if lrv[i] > 0: #right
        #counts
        _ccontra.append(thal_left_counts[i])
        _cipsi.append(thal_right_counts[i])
        #density
        _dcontra.append(thal_density_left[i])
        _dipsi.append(thal_density_right[i])
    elif lrv[i] < 0: #left
        #counts
        _ccontra.append(thal_right_counts[i])
        _cipsi.append(thal_left_counts[i])
        #density
        _dcontra.append(thal_density_right[i])
        _dipsi.append(thal_density_left[i])


_ccontra = np.asarray(_ccontra).T; _dcontra = np.asarray(_dcontra).T
_cipsi = np.asarray(_cipsi).T; _dipsi = np.asarray(_dipsi).T
_dratio = np.asarray([_dcontra[i]/_dipsi[i] for i in range(len(_dcontra))])
_cratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
#make into one
_dist = lrv
 
_inj = np.asarray([inj[i] for i in range(len(inj)) if brains[i] in lr_brains])
_primary_pool = np.asarray([primary_pool[i] for i in range(len(primary_pool)) if brains[i] in lr_brains])

fig = plt.figure()
ax = fig.add_axes([.15,.15,.5,.8])

X = total_cell_counts_per_brain
Y = _cratio.T

ax.scatter(X, Y, color = "k", s = 10)

ax.set_xlabel("Total thalamic cell counts")
ax.set_ylabel("Contra/Ipsi count ratio")
ax.axhline(y=1.2, color = "grey")
plt.savefig(os.path.join(sv_dst, "thal_contra_ipsi_raw_scatter.pdf"), dpi=300)

fig = plt.figure(figsize=(9,5))
ax = fig.add_axes([.15,.15,.5,.8])

#calculate ranks
X = np.argsort(np.sort(total_cell_counts_per_brain))+1
Y = _cratio.T[np.argsort(np.sort(total_cell_counts_per_brain))].argsort()+1 #check if works

results = sm.OLS(Y,sm.add_constant(X)).fit()

m = results.params[1]
r2 = results.rsquared
b = results.params[0]

textstr = "\n".join((
    "slope: {:0.2f}".format(m),
    "$R^2$: {:0.2f}".format(r2)))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

ax.scatter(X, Y, color = "red", s = 10)
ax.plot(m*X + b, color = "k", linestyle = "--")

ytick_spacing = 1; xtick_spacing = 1
ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
ax.set_xlim([0,21])
ax.set_xlabel("Thalamic counts rank order")
ax.set_ylabel("Contra/Ipsi count ratio rank order")
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment="top", bbox=props)
plt.savefig(os.path.join(sv_dst, "thal_contra_ipsi_rank.pdf"), dpi=300)
