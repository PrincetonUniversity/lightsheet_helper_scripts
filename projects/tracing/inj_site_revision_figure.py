#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 13:25:05 2020

@author: wanglab
"""

import os, tifffile, matplotlib.pyplot as plt, numpy as np, matplotlib, pickle, copy
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop

#path to segmented injection sites
# injpth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/injection_sites"
injpth = "/home/wanglab/wang/zahra/tracing_projects/prv/prv_injection_sites"
#path to brain data
brainpth = "/home/wanglab/wang/pisano/tracing_output/retro_4x"
#imports
#path to pickle file
#hsv
# src = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
# data_pth = os.path.join(src, "nc_hsv_maps_contra_pma.p")
# data = pickle.load(open(data_pth, "rb"), encoding = "latin1")
# #set the appropritate variables
# brains = np.array(data["brains"])
# try:
#     frac_of_inj_pool = data["frac_of_inj_pool"]
# except:
#     expr_all_as_frac_of_inj = data["expr_all_as_frac_of_inj"]
#     frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
#                                 for xx in expr_all_as_frac_of_inj])
# ak_pool = data["ak_pool"]

#prv
src = "/jukebox/wang/zahra/tracing_projects/prv"
data_pth = os.path.join(src, "for_tp/prv_maps_contra_pma.p")
data = pickle.load(open(data_pth, "rb"), encoding = "latin1")

#set the appropritate variables
brains = np.array(data["brains"])
frac_of_inj_pool = data["frac_of_inj_pool"]
primary_pool = data["primary_pool"]
ak_pool = data["ak_pool"]

#ideally for the figure we want to show raw data (registered image)
#next to segmentation on top of atlas, next to break by % of injection per lobule
#read atlas
atlpth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
atl = tifffile.imread(atlpth)
atlc = atl[:,460:,:] #crop same way?
#plot segment on top of atlas
atl_cor = np.rot90(np.transpose(atlc, [1,0,2]),axes=(2,1))
    
for brain in os.listdir(injpth):
    
    brain = brain[:-4]
    print(brain)
    try: 
        #find registered injection volume
        regfld = os.path.join(brainpth, brain, "elastix")
        reginjfld = os.path.join(regfld, [xx for xx in os.listdir(regfld) if "555" in xx][0])
        reginjvolpth = os.path.join(reginjfld, "result.tif")
        assert os.path.exists(reginjvolpth)
        #find the injection segmentation of this brain
        seginjvolpth = os.path.join(injpth, brain+".tif")
        assert os.path.exists(seginjvolpth)
        #get all z coordinates from which injection was segmented
        #remember these are sagittal images and cropped in y
        seginjvol = tifffile.imread(seginjvolpth)
        #if not cropped, crop
        if seginjvol.shape == (540, 640, 352):
            seginjvol = seginjvol[:,423:,:]
        nzz,nzy,nzx = np.nonzero(seginjvol)
        
        #now get registered inj volume, crop the same way in y, and maxproject across the range of z-planes in nzz
        reginjvol = tifffile.imread(reginjvolpth)
        # print(reginjvol.shape)
        reginjvolc = reginjvol[:,423:,:]
        #plot
        # plt.imshow(np.max(reginjvolc[np.min(nzz):np.max(nzz)], axis=0)) #sagitta
        #transpose original cropped vol to coronal
        reginjvolc_cor = np.rot90(np.transpose(reginjvolc, [1,0,2]),axes=(2,1))
        #plot
        # plt.imshow(np.max(reginjvolc_cor[np.min(nzy):np.max(nzy)], axis=0)) #coronal
        # plt.imshow(np.max(atl_cor, axis=0), cmap="gist_yarg")
        # plt.imshow(np.max(np.rot90(np.transpose(seginjvol, [1,0,2]),axes=(2,1)), axis=0), cmap=cmap, alpha=0.6)
        
        #get injection fractions
        ind = np.where(brains==brain)[0][0]
        frac_of_inj = frac_of_inj_pool[ind]
        #show
        cmap = copy.copy(plt.cm.Reds)
        cmap.set_over(cmap(1.0))
        cmap.set_under("white")
        vmin = -0.03
        vmax = 0.2
        # plt.imshow([frac_of_inj], cmap=cmap, vmin=vmin, vmax=vmax)
        
        #now put it together in a figure
        #first, separately save out maxprojected registered image you can adjust in imagej
        dst = "/home/wanglab/Desktop/inj_sites"
        tifffile.imsave(os.path.join(dst, brain+"_prv_nc.tif"), 
                np.max(reginjvolc_cor[np.min(nzy):np.max(nzy)], axis=0))
        fig, axes = plt.subplots(ncols = 2, nrows = 1, figsize = (7,5), gridspec_kw = {"wspace":0, "hspace":0,
                                 "width_ratios": [5,0.2]})
        
        #flattened segment on atlas
        ax = axes[0]
        ax.imshow(np.max(atl_cor, axis=0), cmap="gist_yarg")
        ax.imshow(np.max(np.rot90(np.transpose(seginjvol, [1,0,2]),axes=(2,1)), axis=0), 
                  cmap=cmap, alpha=0.8)
        ax.axis("off")
        ax = axes[1]
        #injection site fractions
        show = np.array([frac_of_inj]).T
        cmap = copy.copy(plt.cm.Reds)
        cmap.set_over(cmap(1.0))
        cmap.set_under("white")
        vmin = 0
        vmax = max(show)-0.1
        #colormap
        pc = ax.imshow(show, cmap=cmap, vmin=vmin, vmax=vmax)
        #corresponding color bar
        cb = plt.colorbar(pc, ax=ax, format="%0.1f", orientation="horizontal", shrink=6)#
        cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
        cb.ax.tick_params(labelsize="small")
        cb.ax.set_visible(True) #TP
        ax.set_yticks(np.arange(len(ak_pool)))
        ax.set_xticks([])
        ax.set_yticklabels(ak_pool, fontsize="medium")
        ax.set_xticklabels([])
        ax.tick_params(length=6)
        ax.yaxis.tick_right()
        
        plt.savefig(os.path.join(dst, brain+"_prv_nc.pdf"), dpi=300, bbox_inches="tight") 
        plt.close()
    except:
        print("\nerrors occured\ncheck paths or whether brain doesn't have inj volume\n")                        