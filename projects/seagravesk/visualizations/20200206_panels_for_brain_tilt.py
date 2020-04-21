#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:44:28 2020

@author: wanglab
"""

%matplotlib inline
import tifffile, numpy as np, os, matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

src = "/jukebox/wang/seagravesk/lightsheet/201710_cfos_left_side_only_registration"

brains = ['f37106_mouse1', 'm37109_mouse2', 'm37072_observ', 'f37080_mouse2',
       'm37083_demons', 'm37071_observ', 'f37080_mouse1', 'f37107_demons',
       'f37078_observ', 'f37077_demons', 'm37110_demons', 'f37105_observ', 
       'm37081_observ', 'm37081_demons', 'f37070_observ', 'f37077_observ',
       'm37111_demons', 'm37079_mouse1', 'f37070_demons', 'm37072_demons',
       'f37104_demons', 'm37111_observ', 'f37106_mouse2', 'f37105_demons', 'm37113_mouse1', 'm37071_demons',
       'm37083_observ', 'f37104_observ', 'f37073_mouse1', 'f37107_observ', 'm37113_mouse2', 'm37112_observ',
       'f37073_mouse2', 'm37079_mouse2', 'm37110_observ', 'm37112_demons', 'm37109_mouse1']

flds = [os.path.join(src, xx) for xx in brains]

dst = "/jukebox/wang/zahra/kelly_cell_detection_analysis"

#pick channel and imaging orientation to vis
channel = "790"#"790"
orientation = "dorsal"

#multiply by factor to make images brighter
factor = 20
pdf_pages = PdfPages(os.path.join(dst, "%s_%s_brain_profiles.pdf" % (channel, orientation))) #compiles into multiple pdfs

for i,fld in enumerate(flds[::3]):
    print(fld)
    #brain 1
    arr1_sag = tifffile.imread([os.path.join(fld, xx) for xx in os.listdir(fld) if channel in xx and "resampled" in xx][0])
    arr1_cor = np.transpose(arr1_sag, [1, 2, 0])
    arr1_hor = np.transpose(arr1_sag, [2, 1, 0])
    
    #take max proj of middle planes
    arr1_sag_max = np.max(arr1_sag[250:255], axis = 0)
    arr1_cor_max = np.max(arr1_cor[300:305], axis = 0)
    arr1_hor_max = np.max(arr1_hor[200:205], axis = 0)
    
    #brain 2
    arr2_sag = tifffile.imread([os.path.join(flds[i+1], xx) for xx in os.listdir(flds[i+1]) if channel in xx and "resampled" in xx][0])
    arr2_cor = np.transpose(arr2_sag, [1, 2, 0])
    arr2_hor = np.transpose(arr2_sag, [2, 1, 0])
    
    #take max proj of middle planes
    arr2_sag_max = np.max(arr2_sag[250:255], axis = 0)
    arr2_cor_max = np.max(arr2_cor[300:305], axis = 0)
    arr2_hor_max = np.max(arr2_hor[200:205], axis = 0)
    
    #brain 3
    arr3_sag = tifffile.imread([os.path.join(flds[i+2], xx) for xx in os.listdir(flds[i+2]) if channel in xx and "resampled" in xx][0])
    arr3_cor = np.transpose(arr3_sag, [1, 2, 0])
    arr3_hor = np.transpose(arr3_sag, [2, 1, 0])
    
    #take max proj of middle planes
    arr3_sag_max = np.max(arr3_sag[250:255], axis = 0)
    arr3_cor_max = np.max(arr3_cor[300:305], axis = 0)
    arr3_hor_max = np.max(arr3_hor[200:205], axis = 0)
    
    #plot
    fig = plt.figure(figsize = (15,7), constrained_layout=True)
    gs = fig.add_gridspec(2,6)
    #brain 1
    f_ax1 = fig.add_subplot(gs[1, :2])
    f_ax1.imshow(arr1_cor_max*factor, cmap = "gist_yarg")
    f_ax1.set_title(os.path.basename(fld)+"_"+orientation)
    f_ax1.set_yticks([])
    f_ax1.set_xticks([])
    f_ax2 = fig.add_subplot(gs[0,0])
    f_ax2.imshow(arr1_sag_max*factor, cmap = "gist_yarg")
    f_ax2.set_yticks([])
    f_ax2.set_xticks([])
    f_ax3 = fig.add_subplot(gs[0,1])
    f_ax3.imshow(arr1_hor_max*factor, cmap = "gist_yarg")
    f_ax3.set_yticks([])
    f_ax3.set_xticks([])
    
    #brain 2
    f_ax4 = fig.add_subplot(gs[1, 2:4])
    f_ax4.imshow(arr2_cor_max*factor, cmap = "gist_yarg")
    f_ax4.set_title(os.path.basename(flds[i+1])+"_"+orientation)
    f_ax4.set_yticks([])
    f_ax4.set_xticks([])
    f_ax5 = fig.add_subplot(gs[0,2])
    f_ax5.imshow(arr2_sag_max*factor, cmap = "gist_yarg")
    f_ax5.set_yticks([])
    f_ax5.set_xticks([])
    f_ax6 = fig.add_subplot(gs[0,3])
    f_ax6.imshow(arr2_hor_max*factor, cmap = "gist_yarg")
    f_ax6.set_yticks([])
    f_ax6.set_xticks([])
    
    #brain 3
    f_ax7 = fig.add_subplot(gs[1, 4:6])
    f_ax7.imshow(arr3_cor_max*factor, cmap = "gist_yarg")
    f_ax7.set_title(os.path.basename(flds[i+2])+"_"+orientation)
    f_ax7.set_yticks([])
    f_ax7.set_xticks([])
    f_ax8 = fig.add_subplot(gs[0,4])
    f_ax8.imshow(arr3_sag_max*factor, cmap = "gist_yarg")
    f_ax8.set_yticks([])
    f_ax8.set_xticks([])
    f_ax9 = fig.add_subplot(gs[0,5])
    f_ax9.imshow(arr3_hor_max*factor, cmap = "gist_yarg")
    f_ax9.set_yticks([])
    f_ax9.set_xticks([])
    
    
    #done with the page
    pdf_pages.savefig(dpi = 300, bbox_inches = "tight") 
    plt.close()

#write PDF document contains all brains
pdf_pages.close()