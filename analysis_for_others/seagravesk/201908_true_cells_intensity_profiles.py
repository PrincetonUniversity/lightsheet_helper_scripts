#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:03:34 2019

@author: wanglab
"""

import os, numpy as np, tifffile, matplotlib.pyplot as plt, cv2
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages
from tools.conv_net.utils.io import read_roi_zip
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop

"""
One solution that was suggested to me was to check every putative cell center for whether it is a peak in 
all dimensions (i.e. z, x, and y). The idea is that a cell should be roughly spherical in shape, whereas edges will not be.

For example, should we take the three orthogonal vectors that pass through the voxel labeled as the cell center 
and do peak detection on those vectors?

Number of cells in Volume 1: 126
Number of cells in Volume 2: 165
"""

src = "/jukebox/wang/seagravesk/lightsheet/man_label_cells_Kelly/In prog/complete"
raw = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos"
w = 10
dst = "/home/wanglab/Desktop/real_cell_profiles.pdf"
pdf_pages = PdfPages(dst) #compiles into multiple pdfs

#volume 1
impths = [os.path.join(src, "171209_f37104_demons/171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13_plane500to519_volume1.tif"),
         os.path.join(src, "171210_m37079_mouse2/171210_m37079_mouse2_20171014_790_015na_1hfsds_z5um_1000msec_17-25-29_planes280to299_volume1.tif")]

fullszimpth = [os.path.join(raw, "171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13"),
               os.path.join(raw, "171210_m37079_mouse2_20171014_790_015na_1hfsds_z5um_1000msec_17-25-29")]

roi_pths = [os.path.join(src, "171209_f37104_demons/171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13_plane500to519_volume1_v11.RoiSet.zip"),
           os.path.join(src, "171210_m37079_mouse2/171210_m37079_mouse2_20171014_790_015na_1hfsds_z5um_1000msec_17-25-29_planes280to299_volume1_v9.RoiSet.zip")]

fullszbnds = [(500, 590, 937), (280, 893, 934)] #from kelly notes

#make dict of cell coordinates (IN FIJI TERMS) and profile(s) to save
pnts_dct = {}

for i, impth in enumerate(impths):
    
    vol = tifffile.imread(impth)
    if i == 0:
        fullszstkpth = [os.path.join(fullszimpth[i], xx) for xx in os.listdir(fullszimpth[i]) if "Z049" in xx or "Z050" #custom
                        in xx or "Z051" in xx or "Z052" in xx]; fullszstkpth.sort() #40 planes
    else:
        fullszstkpth = [os.path.join(fullszimpth[i], xx) for xx in os.listdir(fullszimpth[i]) if "Z027" in xx or "Z028" #custom
                        in xx or "Z029" in xx or "Z030" in xx]; fullszstkpth.sort() #40 planes
        
    fullszstk = np.asarray([tifffile.imread(pth) for pth in fullszstkpth])
    #note python 3 change
    zyx_rois = np.asarray([[int(yy) for yy in xx[0].replace(".roi", "").split("-")] for xx in
                            read_roi_zip(roi_pths[i], include_roi_name=True)])
    
    pnts = np.asarray([[xx[0]-1, xx[1], xx[2]] for xx in zyx_rois])
    bnds = fullszbnds[i]
    fullsz_pnts = np.asarray([[xx[0]+9, xx[1]+bnds[1], xx[2]+bnds[2]] for xx in zyx_rois])
    
    #save all features of all points in volume
    big_dct = {}
    for i, pntn in enumerate(pnts):
        
        #save features to dict
        dct = {}        
        print(pntn)
        
        fszpnt = fullsz_pnts[i]
        seg = fullszstk[fszpnt[0], fszpnt[1]-w:fszpnt[1]+w, fszpnt[2]-w:fszpnt[2]+w]
        #have done these checks, obv doesnt work on edge volumes
#        ann_seg = vol[pntn[0], pntn[1]-w:pntn[1]+w, pntn[2]-w:pntn[2]+w]
#        assert np.sum(seg-ann_seg)==0 #make sure its the same image
        
        #map cell map
        cell = np.zeros_like(seg)
        cell[int(cell.shape[0]/2), int(cell.shape[1]/2)] = 100        
        kernel = np.ones((1,1), np.uint8) 
        cell_dilation = cv2.dilate(cell, kernel, iterations=1) 
    
        #saving
        fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (6,7), gridspec_kw = {"wspace":0, "hspace":0,
                                 "height_ratios": [2,1]})
        
        intensity = vol[pntn[0], pntn[1], pntn[2]]
        axes[0].imshow(seg*20, "gist_yarg")
        axes[0].imshow(cell_dilation, cmap, alpha = 0.6)
        axes[0].set_ylabel("Location of real cell")
        axes[0].set_title("Volume #{}\nCell coordinate (zyx): {}\nIntensity = {}".format(i+1, tuple(pntn), intensity))
        axes[0].set_xticks([], [])
        axes[0].set_yticks([], [])
        
        yprofile = fullszstk[fszpnt[0], fszpnt[1]-w:fszpnt[1]+w+1, fszpnt[2]] #cell center indexed at 10
        xprofile = fullszstk[fszpnt[0], fszpnt[1], fszpnt[2]-w:fszpnt[2]+w+1]
        zprofile = fullszstk[fszpnt[0]-w:fszpnt[0]+w+1, fszpnt[1], fszpnt[2]]
            
        axes[1].plot(zprofile, label = "Z")
        axes[1].plot(yprofile, label = "Y")
        axes[1].plot(xprofile, label = "X")
        axes[1].set_ylabel("Pixel intensity")
        axes[1].set_xlabel("Normalized distance")
        axes[1].set_xticks(np.arange(0, len(xprofile), 2))
        axes[1].set_xticklabels([-1. , -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8], fontsize = "x-small")
        axes[1].legend()
    
        #done with the page
        pdf_pages.savefig(dpi = 300, bbox_inches = 'tight') 
        plt.close()
        
        print("saving...")
        #save profiles to dct
        dct["yprofile"] = yprofile
        dct["xprofile"] = xprofile
        dct["zprofile"] = zprofile
        dct["cell_center"] = tuple(pntn)
        dct["intensity"] = intensity
        dct["img_segment"] = seg
        
        big_dct[tuple(pntn)] = dct

    #save points per volume
    pnts_dct[os.path.basename(os.path.dirname(impth))] = big_dct 
    

#write PDF document contains all points
pdf_pages.close()

print("\nexporting to pickle...\n")
#save data to pickle
import pickle
sv = "/jukebox/wang/zahra/kelly_cell_detection_analysis/real_cells.p"

with open(sv, "wb") as fp:
    pickle.dump(pnts_dct, fp, protocol=pickle.HIGHEST_PROTOCOL)

data = pickle.load(open(sv, "rb"), encoding = "latin1")