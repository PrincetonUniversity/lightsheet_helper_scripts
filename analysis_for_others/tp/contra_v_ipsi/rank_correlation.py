#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 17:21:19 2019

@author: wanglab
"""

import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, matplotlib as mpl

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


fig = plt.figure()
ax = fig.add_axes([.15,.15,.5,.8])

X = np.sort(total_cell_counts_per_brain)
Y = _cratio.T[np.argsort(total_cell_counts_per_brain)]

ax.scatter(X, Y, color = "k", s = 10)

ax.set_xlabel("Total cell counts")
ax.set_ylabel("Contra/Ipsi count ratio")
ax.axhline(y=1.2, color = "grey", linestyle = "--")
plt.savefig(os.path.join(sv_dst, "contra_ipsi.pdf"), dpi=300)