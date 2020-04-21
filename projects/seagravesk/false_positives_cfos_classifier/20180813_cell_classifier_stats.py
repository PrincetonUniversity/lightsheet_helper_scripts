#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:12:37 2019

@author: wanglab
"""

import pickle, numpy as np, matplotlib.pyplot as plt, os

src = "/jukebox/wang/zahra/kelly_cell_detection_analysis"
ecells = os.path.join(src, "edge_cells.p")

rcells = os.path.join(src, "real_cells.p")

rcells = pickle.load(open(rcells, "rb"), encoding = "latin1")
ecells = pickle.load(open(ecells, "rb"), encoding = "latin1")

#change to see distance in pixels around cell center
cutoff = 3
#look at real cells first
#do for all x,y,z
yprof = []
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        yprof.append(cdct[j]["yprofile"])
#make np array
yprof = np.asarray(yprof)
#normalize
norm_yprof = []
for yy in yprof:    
    norm_yprof.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))    
norm_yprof = np.asarray(norm_yprof)

xprof = []
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        xprof.append(cdct[j]["xprofile"])
#make np array
xprof = np.asarray(xprof)
#normalize
norm_xprof = []
for yy in xprof:    
    norm_xprof.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_xprof = np.asarray(norm_xprof)

zprof = []
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        zprof.append(cdct[j]["zprofile"])
#make np array
zprof = np.asarray(zprof)
#normalize
norm_zprof = []
for yy in zprof:    
    norm_zprof.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_zprof = np.asarray(norm_zprof)

#%%
#do the same for edge cells

#do for all x,y,z
yprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    if type(v) is dict:
        yprof_e.append(cdct["yprofile"])
#make np array
yprof_e = np.asarray(yprof_e)
#normalize
norm_yprof_e = []
for yy in yprof_e:    
    norm_yprof_e.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))    
norm_yprof_e = np.asarray(norm_yprof_e)

xprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    if type(v) is dict:
        xprof_e.append(cdct["xprofile"])
#make np array
xprof_e = np.asarray(xprof_e)
#normalize
norm_xprof_e = []
for yy in xprof_e:    
    norm_xprof_e.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_xprof_e = np.asarray(norm_xprof_e)

zprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    if type(v) is dict:
        zprof_e.append(cdct["zprofile"])
#make np array
zprof_e = np.asarray(zprof_e)
#normalize
norm_zprof_e = []
for yy in zprof_e:    
    norm_zprof_e.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_zprof_e = np.asarray(norm_zprof_e)

#%%
#make into sensible numpy arrays to take mean
norm_z = np.ones((norm_zprof.shape[0], 21))*float("nan")
for i in range(len(norm_zprof)):
    norm_z[i, :len(norm_zprof[i])] = norm_zprof[i]
#mean ignoring nans and last 2 values
norm_z_mean = np.nanmean(norm_z[:,:19], axis = 0)    
norm_y = np.ones((norm_yprof.shape[0], 21))*float("nan")
for i in range(len(norm_yprof)):
    norm_y[i, :len(norm_yprof[i])] = norm_yprof[i]
norm_y_mean = np.nanmean(norm_y, axis = 0)        
norm_x = np.ones((norm_xprof.shape[0], 21))*float("nan")
for i in range(len(norm_xprof)):
    norm_x[i, :len(norm_xprof[i])] = norm_xprof[i]
norm_x_mean = np.nanmean(norm_x, axis = 0)        
#%%
#plotting and saving
# plot of all the overlays + mean 
fig, axes = plt.subplots(ncols = 2, nrows = 3, figsize = (6,7), sharex = True, sharey = True, gridspec_kw = {"wspace":0, "hspace":0,
                                     "height_ratios": [1,1,1]})
lw = 0.5
alpha = 0.5
#plot all normalized
for i in range(len(norm_zprof)):
    if i == 0:
        axes[0,0].plot(norm_zprof[i], linewidth = lw, color = "gray", linestyle = "dashed", alpha = alpha)
        axes[0,0].plot(norm_z_mean, linewidth = lw*10, color = "k", label = "Z (mean)")
    else:
        axes[0,0].plot(norm_zprof[i], linewidth = lw, color = "gray", linestyle = "dashed", alpha = alpha)
axes[0,0].set_ylabel("Normalized pixel intensity", fontsize = "x-small")

axes[0,0].legend()
axes[0,0].set_title("Human annotated cells\n(n={})".format(len(norm_xprof)))

for i in range(len(norm_yprof)):
    if i == 0:
        axes[1,0].plot(norm_yprof[i], linewidth = lw, color = "lightcoral", linestyle = "dashed", alpha = alpha)
        axes[1,0].plot(norm_y_mean, linewidth = lw*10, color = "red", label = "Y (mean)")
    else:
        axes[1,0].plot(norm_yprof[i], linewidth = lw, color = "lightcoral", linestyle = "dashed", alpha = alpha)
axes[1,0].set_ylabel("Normalized pixel intensity", fontsize = "x-small")

axes[1,0].legend()
for i in range(len(norm_xprof)):
    if i == 0:
        axes[2,0].plot(norm_xprof[i], linewidth = lw, color = "lightsteelblue", linestyle = "dashed", alpha = alpha)
        axes[2,0].plot(norm_x_mean, linewidth = lw*10, color = "blue", label = "X (mean)")
    else:
        axes[2,0].plot(norm_xprof[i], linewidth = lw, color = "lightsteelblue", linestyle = "dashed", alpha = alpha)
axes[2,0].set_xticks(np.arange(0, len(norm_xprof[0])+1, 2))
axes[2,0].set_xlabel("Distance (pixels)")
axes[2,0].set_ylabel("Normalized pixel intensity", fontsize = "x-small")

axes[2,0].legend()

#plot all normalized
for i in range(len(norm_zprof_e)):
    if i == 0:
        axes[0,1].plot(np.mean(norm_zprof_e, axis = 0), linewidth = lw*10, color = "k")
    else:
        axes[0,1].plot(norm_zprof_e[i], linewidth = lw, color = "gray", linestyle = "dashed", alpha = alpha)

axes[0,1].set_title("'Edge' cells\n(n={})".format(len(norm_xprof_e)))

for i in range(len(norm_yprof_e)):
    if i == 0:
        axes[1,1].plot(np.mean(norm_yprof_e,axis = 0), linewidth = lw*10, color = "red")
    else:
        axes[1,1].plot(norm_yprof_e[i], linewidth = lw, color = "lightcoral", linestyle = "dashed", alpha = alpha)

for i in range(len(norm_xprof_e)):
    if i == 0:
        axes[2,1].plot(np.mean(norm_xprof_e,axis = 0), linewidth = lw*10, color = "blue")
    else:
        axes[2,1].plot(norm_xprof_e[i], linewidth = lw, color = "lightsteelblue", linestyle = "dashed", alpha = alpha)
    
axes[2,1].set_xticks(np.arange(0, len(norm_xprof_e[0])+1, 2))
axes[2,1].set_xlabel("Distance (pixels)")

plt.savefig(os.path.join(src, "overlays_real_edge_cells_w_mean.pdf"), dpi = 300, bbox_inches = "tight")

#%%
#big plot of all the overlays
fig, axes = plt.subplots(ncols = 2, nrows = 3, figsize = (6,7), sharex = True, sharey = True, gridspec_kw = {"wspace":0, "hspace":0,
                                     "height_ratios": [1,1,1]})
lw = 0.5
#plot all normalized
for i in range(len(norm_zprof)):
    if i == 0:
        axes[0,0].plot(norm_zprof[i], linewidth = lw, color = "k", label = "Z")
    else:
        axes[0,0].plot(norm_zprof[i], linewidth = lw, color = "k")
axes[0,0].set_ylabel("Normalized pixel intensity", fontsize = "x-small")

axes[0,0].legend()
axes[0,0].set_title("Human annotated cells\n(n={})".format(len(norm_xprof)))

for i in range(len(norm_yprof)):
    if i == 0:
        axes[1,0].plot(norm_yprof[i], linewidth = lw, color = "red", label = "Y")
    else:
        axes[1,0].plot(norm_yprof[i], linewidth = lw, color = "red")
axes[1,0].set_ylabel("Normalized pixel intensity", fontsize = "x-small")

axes[1,0].legend()
for i in range(len(norm_xprof)):
    if i == 0:
        axes[2,0].plot(norm_xprof[i], linewidth = lw, color = "blue", label = "X")
    else:
        axes[2,0].plot(norm_xprof[i], linewidth = lw, color = "blue")
axes[2,0].set_xticks(np.arange(0, len(norm_xprof[0])+1, 2))
axes[2,0].set_xlabel("Distance (pixels)")
axes[2,0].set_ylabel("Normalized pixel intensity", fontsize = "x-small")

axes[2,0].legend()

#plot all normalized
for i in range(len(norm_zprof_e)):
    axes[0,1].plot(norm_zprof_e[i], linewidth = lw, color = "k")

axes[0,1].set_title("'Edge' cells\n(n={})".format(len(norm_xprof_e)))

for i in range(len(norm_yprof_e)):
    axes[1,1].plot(norm_yprof_e[i], linewidth = lw, color = "red")

for i in range(len(norm_xprof_e)):
    axes[2,1].plot(norm_xprof_e[i], linewidth = lw, color = "blue")
axes[2,1].set_xticks(np.arange(0, len(norm_xprof_e[0])+1, 2))
axes[2,1].set_xlabel("Distance (pixels)")

plt.savefig(os.path.join(src, "overlays_real_edge_cells.pdf"), dpi = 300, bbox_inches = "tight")