#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:59:53 2019

@author: wanglab
"""

import os, tifffile, numpy as np, matplotlib.pyplot as plt
import matplotlib.colors

from tools.utils.io import load_kwargs

def mask_atlas_w_coordinates(atl_pth, atl_id):
    """ returns z,y,x coordinates of structure requested """
    
    atl = tifffile.imread(atl_pth)
    #orient horizontally
    atl_hor = np.transpose(atl, [2, 1, 0])
    #find zyx coordinates of structure
    zplns,y,x = np.where(atl_hor == atl_id)
    #make mask
    atl_hor[atl_hor != atl_id] = 0
    
    return zplns,y,x, atl_hor
    
if __name__ == "__main__":
    
    src = "/jukebox/wang/pisano/tracing_output/antero_4x"
    dst = "/home/wanglab/Desktop/pons_mask_figs"
    atl_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif"
    
    #animals = ["an01", "an02", "an03", "an04", "an05", "an06", "an07", "an09", "an10", "an12", "an13", "an14", "an15", "an16", "an17"]
    animals = ["20170410_tp_bl6_lob6a_ml_repro_01", "20160823_tp_bl6_cri_500r_02", "20180417_jg59_bl6_cri_03",
            "20170207_db_bl6_crii_1300r_02", "20160622_db_bl6_unk_01", "20161205_tp_bl6_sim_750r_03",
            "20180410_jg51_bl6_lob6b_04", "20170419_db_bl6_cri_rpv_53hr", "20170116_tp_bl6_lob6b_lpv_07",
            "20170411_db_bl6_crii_mid_53hr", "20160822_tp_bl6_crii_1500r_06", "20160920_tp_bl6_lob7_500r_03",
            "20170207_db_bl6_crii_rpv_01", "20161205_tp_bl6_sim_250r_02", "20161207_db_bl6_lob6a_500r_53hr",
            "20170130_tp_bl6_sim_rlat_05", "20170115_tp_bl6_lob6b_500r_05", "20170419_db_bl6_cri_mid_53hr",
            "20161207_db_bl6_lob6a_850r_53hr", "20160622_db_bl6_crii_52hr_01", "20161207_db_bl6_lob6a_50rml_53d5hr",
            "20161205_tp_bl6_lob45_1000r_01", "20160801_db_l7_cri_01_mid_64hr"]
    
    brains = [os.path.join(src, xx) for xx in animals]
    
    #make binary colarmap
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "lime"])
    
    #setting for masking regions
    maxip_scale = 4 #aka 60 um sections
    atl_id = 771 #mask pons
    zplns,y,x, atl_hor = mask_atlas_w_coordinates(atl_pth, atl_id)
    maxip_stop = max(zplns)-50
    maxip_start = min(zplns)+50
    yrange = (min(y)-10, max(y)+10) #10 is the padding around the edge of the structure
    xrange = (min(x)-10, max(x)+10) 
    
    for brain in brains:
        
        kwargs = load_kwargs(brain)
        cellvol = [xx for xx in kwargs["volumes"] if xx.ch_type == "cellch"][0]
        vol = tifffile.imread(cellvol.ch_to_reg_to_atlas)
        vol_hor = np.transpose(vol, [2, 1, 0])
        
        #init figure
        f = np.sqrt((maxip_stop - maxip_start)/maxip_scale).astype(int) #has to be a square 
        ncols, nrows = f, f
        fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = (15,10), sharex = True, gridspec_kw = {"wspace":0, "hspace":0})
        slcs = np.arange(maxip_start, maxip_stop, maxip_scale)
        k = 0 #init slice range
        for i in range(ncols):
            for j in range(nrows):
                #crop after max projection
                axes[i,j].imshow(np.max(vol_hor[slcs[k]:slcs[k]+maxip_scale]*8,
                    axis=0)[yrange[0]:yrange[1], xrange[0]:xrange[1]], cmap="Greys")
                axes[i,j].imshow(np.max(atl_hor[slcs[k]:slcs[k]+maxip_scale], 
                    axis=0)[yrange[0]:yrange[1], xrange[0]:xrange[1]], alpha = 0.2, cmap=cmap)
                axes[i,j].axis("off")
                k += 1        
        
        #done with the page
        plt.savefig(os.path.join(dst, os.path.basename(brain)+"_h129_pons_z%d-%d_zstep%d.png" % (maxip_start, 
                                          maxip_stop, maxip_scale)), bbox_inches = "tight")             
        
        plt.close()