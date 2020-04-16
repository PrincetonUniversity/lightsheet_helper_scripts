#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:46:44 2020

@author: wanglab
"""


import os, subprocess as sp, tifffile, numpy as np, matplotlib.pyplot as plt, matplotlib as mpl, pandas as pd
from collections import Counter
from tools.imageprocessing.orientation import fix_orientation

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

if __name__ == "__main__":

    src = "/jukebox/wang/pisano/tracing_output/retro_4x"
    
    #thalamic
    brains = ["20170410_tp_bl6_lob6a_ml_repro_01", "20160823_tp_bl6_cri_500r_02", "20180417_jg59_bl6_cri_03",
        "20170207_db_bl6_crii_1300r_02", "20160622_db_bl6_unk_01", "20161205_tp_bl6_sim_750r_03",
        "20180410_jg51_bl6_lob6b_04", "20170419_db_bl6_cri_rpv_53hr", "20170116_tp_bl6_lob6b_lpv_07",
        "20170411_db_bl6_crii_mid_53hr", "20160822_tp_bl6_crii_1500r_06", "20160920_tp_bl6_lob7_500r_03",
        "20170207_db_bl6_crii_rpv_01", "20161205_tp_bl6_sim_250r_02", "20161207_db_bl6_lob6a_500r_53hr",
        "20170130_tp_bl6_sim_rlat_05", "20170115_tp_bl6_lob6b_500r_05", "20170419_db_bl6_cri_mid_53hr",
        "20161207_db_bl6_lob6a_850r_53hr", "20160622_db_bl6_crii_52hr_01", "20161207_db_bl6_lob6a_50rml_53d5hr",
        "20161205_tp_bl6_lob45_1000r_01", "20160801_db_l7_cri_01_mid_64hr"]    
    
    #neocortical 
    # brains = ["20180409_jg46_bl6_lob6a_04", "20180608_jg75",
    #    "20170204_tp_bl6_cri_1750r_03", "20180608_jg72",
    #    "20180416_jg56_bl6_lob8_04", "20170116_tp_bl6_lob45_ml_11",
    #    "20180417_jg60_bl6_cri_04", "20180410_jg52_bl6_lob7_05",
    #    "20170116_tp_bl6_lob7_1000r_10", "20180409_jg44_bl6_lob6a_02",
    #    "20180410_jg49_bl6_lob45_02", "20180410_jg48_bl6_lob6a_01",
    #    "20180612_jg80", "20180608_jg71", "20170212_tp_bl6_crii_1000r_02",
    #    "20170115_tp_bl6_lob6a_rpv_03", "20170212_tp_bl6_crii_2000r_03",
    #    "20180417_jg58_bl6_sim_02", "20170130_tp_bl6_sim_1750r_03",
    #    "20170115_tp_bl6_lob6b_ml_04", "20180410_jg50_bl6_lob6b_03",
    #    "20170115_tp_bl6_lob6a_1000r_02", "20170116_tp_bl6_lob45_500r_12",
    #    "20180612_jg77", "20180612_jg76", "20180416_jg55_bl6_lob8_03",
    #    "20170115_tp_bl6_lob6a_500r_01", "20170130_tp_bl6_sim_rpv_01",
    #    "20170204_tp_bl6_cri_1000r_02", "20170212_tp_bl6_crii_250r_01",
    #    "20180417_jg61_bl6_crii_05", "20170116_tp_bl6_lob7_ml_08",
    #    "20180409_jg47_bl6_lob6a_05"]
    
    #run
    inputlist = [os.path.join(src, xx) for xx in brains]
    
    dct = {"inputlist": inputlist,
      "channel": "01",
      "channel_type": "injch",
      "filter_kernel": (3,3,3), 
      "threshold": 4, 
      "num_sites_to_keep": 1,
      "injectionscale": 45000, 
      "imagescale": 2,
      "reorientation": ("2","0","1"),
      "crop": "[:,:,:]",
      "dst": "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites",
      "save_individual": False, 
      "save_tif": True,
      "colormap": "plasma", 
      "atlas": "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif",
      "annotation": "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif",
      "id_table": "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"
    }
    
    #inspect injection sites for the brains i currently have
    imgs = [os.path.join(dct["dst"], xx+".tif.tif") for xx in brains]; imgs.sort()
    #the y-axis cutoff for visualization
    sites = np.array([fix_orientation(tifffile.imread(xx), dct["reorientation"]) for xx in imgs]) 
    atl = fix_orientation(tifffile.imread(dct["atlas"])[:, 450:, :], dct["reorientation"]) #cut at 423 for HSV??
    #make counter
    nonzeros = [list(zip(*np.nonzero(site))) for site in sites] #<-for pooled image
    
    #condense nonzero pixels
    nzs = [str(x) for xx in nonzeros for x in xx] #this list has duplicates if two brains had the same voxel w label
    c = Counter(nzs)
    array = np.zeros(atl.shape)
    print("Collecting nonzero pixels for pooled image...")
    tick = 0
    #generating pooled array where voxel value = total number of brains with that voxel as positive
    for k,v in c.items():
        k = [int(xx) for xx in k.replace("(","").replace(")","").split(",")]
        array[k[0], k[1], k[2]] = int(v)
        tick+=1
        if tick % 50000 == 0: print("   {}".format(tick))
        
    my_cmap = eval("plt.cm.{}(np.arange(plt.cm.Reds.N))".format("Reds"))
    my_cmap[:1,:4] = 0.0  
    my_cmap = mpl.colors.ListedColormap(my_cmap)
    my_cmap.set_under("w")
    plt.figure()
    plt.imshow(np.max(atl, axis=0), cmap="gray_r")
    plt.imshow(np.max(array, axis=0), alpha=0.90, cmap=my_cmap); plt.colorbar(); plt.axis("off")
    plt.tick_params(length=6)

    plt.savefig("/home/wanglab/Desktop/hsv_thal_inj_heatmap_cb_inverted_reds.pdf", dpi = 300, transparent = True);
    plt.close()