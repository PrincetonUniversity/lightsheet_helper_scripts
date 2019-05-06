#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 15:39:12 2019

@author: wanglab
"""

import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, os
import SimpleITK as sitk
from skimage.external import tifffile
import numpy as np, scipy
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.registration.allen_structure_json_to_pandas import annotation_location_to_structure
from tools.analysis.network_analysis import make_structure_objects
import matplotlib
tab20cmap = [plt.cm.tab20(xx) for xx in range(20)]
tab20cmap_nogray = tab20cmap[:14] + tab20cmap[16:]
import matplotlib as mpl
tab20cmap_nogray = matplotlib.colors.ListedColormap(tab20cmap_nogray, name = "tab20cmap_nogray")
import matplotlib.patches as mpatches

def consolidate_parents_structures_cfos(id_table, ann, namelist, verbose=False, structures=False):
    """Function that generates evenly spaced pixels values based on annotation parents

    Removes 0 from list

    Inputs:
        id_table=path to excel file generated from scripts above
        ann = allen annoation file
        namelist=list of structues names, typically parent structures*********************

    Returns:
        -----------
        nann = new array of bitdepth
        list of value+name combinations
    """
    if type(ann) == str: ann = sitk.GetArrayFromImage(sitk.ReadImage(ann))


    #remove duplicates and null and root
    namelist = list(set(namelist))
    namelist = [xx for xx in namelist if xx != "null" and xx != "root"]
    namelist.sort()

    #make structures to find parents
    if not structures:
        structures = make_structure_objects(id_table)

    #setup
    nann = np.zeros(ann.shape).astype("uint8")
    cmap = [xx for xx in np.linspace(1,255, num=len(namelist))]

    #populate
    for i in range(len(namelist)):
        try:
            nm=namelist[i]
            s = [xx for xx in structures if xx.name==nm][0]
            if verbose: print ("{}, {} of {}, value {}".format(nm, i, len(namelist)-1, cmap[i]))
            nann[np.where(ann==int(s.idnum))] = cmap[i]
            for ii in s.progeny:
                if ii[3] != "null": nann[np.where(ann==int(ii[3]))] = cmap[i]
        except Exception,e:
            print nm, e
    #sitk.Show(sitk.GetImageFromArray(nann))
    #change nann to have NAN where zeros
    nann = nann.astype("float")
    nann[nann == 0] = "nan"

    return nann, zip(cmap[:], namelist)

def make_2D_overlay_of_heatmaps(src, atl_pth, ann_pth, allen_id_table, save_dst, condition,
                                positive = True, negative = True):
    
    """ custom script for cfos visualisation """
    print(save_dst)

    #half brain        
    #read variables
    vol = tifffile.imread(src)
    atl = tifffile.imread(atl_pth)
    ann = tifffile.imread(ann_pth)
    
    pvol = vol[:,:,:,1]
    nvol = vol[:,:,:,0]
    atl = np.rot90(np.transpose(atl, [1, 0, 2]), axes = (2,1)) #sagittal to coronal
    ann = np.rot90(np.transpose(ann, [1, 0, 2]), axes = (2,1)) #sagittal to coronal
    
    print(atl.shape)
    print(ann.shape)
    print(nvol.shape)
    
    #threshold values
    pvol[pvol!=0.0] = 1.0
    nvol[nvol!=0.0] = 1.0
    
    #loop only positives
    zstep = 50
    colorbar_cutoff = 65#20,60 #this is related to zstep size...(it"s like  apercetnage...)
    rngs = range(50, 500, zstep)
    
    #positive
    if positive:
        print("positive")
        #make fig
        #pvals
        plt.style.use('dark_background')
        fig, axs = plt.subplots(3,3, sharex='col', sharey='row',
                    gridspec_kw={'hspace': 0, 'wspace': 0}, facecolor="black",
                    figsize = (6, 4))                    
        #modify colormap
        my_cmap = plt.cm.viridis(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under("w")       
        #plot
        i = 0
        for row in axs:
            for col in row:
                try:
                    img = col.imshow(np.max(atl[rngs[i]:rngs[i+1]], axis=0), cmap="gray", alpha=1)
                    img = col.imshow(np.sum(pvol[rngs[i]:rngs[i+1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
                    col.axis("off")
                    plt.tight_layout()
                except:
                     col.imshow(np.zeros_like(np.max(atl[rngs[0]:rngs[1]], axis=0)), cmap = "binary_r")
                     col.axis("off")
                     plt.tight_layout()    
                i += 1        
        #custom colorbar
        cbar = fig.colorbar(img, ax=axs.ravel().tolist(), fraction=0.015, pad=0.01)
        cbar.set_label('# of planes represented',size=7)
        # access to cbar tick labels:
        cbar.ax.tick_params(labelsize=5) 
#        cbar.ax.set_yticklabels(['5', '10', '15', '20'])  # vertically oriented colorbar
        
        plt.savefig(os.path.join(save_dst, "{}_positively_correlated_voxels_AP_dimension_half_brain.pdf".format(condition)), bbox_inches = 'tight', dpi=300, 
                    facecolor=fig.get_facecolor(), edgecolor='none', transparent = True)
        plt.close()

    #negative
    if negative:        
        print("negative")        
        #make fig
        #pvals
        plt.style.use('dark_background')
        fig, axs = plt.subplots(3,3, sharex='col', sharey='row',
                    gridspec_kw={'hspace': 0, 'wspace': 0}, facecolor="black",
                    figsize = (6, 4))                             
        #modify colormap
        my_cmap = plt.cm.plasma(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under("w")        
        #plot
        i = 0
        for row in axs:
            for col in row:
                try:
                    img = col.imshow(np.max(atl[rngs[i]:rngs[i+1]], axis=0), cmap="gray", alpha=1)
                    img = col.imshow(np.sum(nvol[rngs[i]:rngs[i+1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
                    col.axis("off")
                    plt.tight_layout()
                except:
                     col.imshow(np.zeros_like(np.max(atl[rngs[0]:rngs[1]], axis=0)), cmap = "binary_r")
                     col.axis("off")  
                     plt.tight_layout()
                i += 1
        
        cbar = fig.colorbar(img, ax=axs.ravel().tolist(), fraction=0.015, pad=0.01)
        cbar.set_label('# of planes represented',size=7)
        # access to cbar tick labels:
        cbar.ax.tick_params(labelsize=5) 
#        cbar.ax.set_yticklabels(['5', '10', '15', '20'])  # vertically oriented colorbar
        
        plt.savefig(os.path.join(save_dst, "{}_negatively_correlated_voxels_AP_dimension_half_brain.pdf".format(condition)), bbox_inches = 'tight', dpi=300, 
                    facecolor=fig.get_facecolor(), edgecolor='none', transparent = True)
        plt.close()
       
    return 

#%%
if __name__ == "__main__":

    #set destination of p value map you want to analyze
    ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_60um_erosion.tif"
    atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
    allen_id_table = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
    dst = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/p_value_maps"    
    save_dst = os.path.join(dst, "composite_stacks")
    if not os.path.exists(save_dst): os.mkdir(save_dst)
    
    #REMEMBER THAN THE FIRST GROUP IS THE CONTROL = RED, SECOND GROUP IS STIMULATION = BLUE
    #run
    #condition 1
    src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/p_value_maps/pvalues_CNO_control_no_reversal_vs_CNO_control_reversal.tif"            
    #make parent list and 2d overlays
    make_2D_overlay_of_heatmaps(src, atl_pth, ann_pth, allen_id_table, save_dst, 
                                condition = "CNO_control_no_reversal_vs_CNO_control_reversal", positive = True, negative = True)
    
    #condition 2
    src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/p_value_maps/pvalues_homecage_control_vs_CNO_control_no_reversal.tif"
    make_2D_overlay_of_heatmaps(src, atl_pth, ann_pth, allen_id_table, save_dst, 
                                condition = "homecage_control_vs_CNO_control_no_reversal", positive = True, negative = True)
    
    #condition 3
    src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/p_value_maps/pvalues_CNO_control_reversal_vs_DREADDs.tif"
    make_2D_overlay_of_heatmaps(src, atl_pth, ann_pth, allen_id_table, save_dst, 
                                condition = "CNO_control_reversal_vs_DREADDs", positive = True, negative = True)
    
    #condition 4
    src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/p_value_maps/pvalues_homecage_control_vs_DREADDs.tif"
    make_2D_overlay_of_heatmaps(src, atl_pth, ann_pth, allen_id_table, save_dst, 
                                condition = "homecage_control_vs_DREADDs", positive = True, negative = True)
    
    
    