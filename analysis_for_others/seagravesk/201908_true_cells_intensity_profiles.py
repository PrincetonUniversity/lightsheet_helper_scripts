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

Unfortunately, I don’t really know how to do this, which is where I need help. :-)
I’m sure that this is a common problem in image processing / image analysis, so it is likely just a matter 
of looking into what other people have done.

For example, should we take the three orthogonal vectors that pass through the voxel labeled as the cell center 
and do peak detection on those vectors?

Number of cells in Volume 1: 126

Number of cells in Volume 2: 165

"""

src = "/jukebox/wang/seagravesk/lightsheet/man_label_cells_Kelly/In prog/complete"

w = 10
dst = "/home/wanglab/Desktop/real_cell_profiles.pdf"
pdf_pages = PdfPages(dst) #compiles into multiple pdfs

#volume 1
impths = [os.path.join(src, "171209_f37104_demons/171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13_plane500to519_volume1.tif"),
         os.path.join(src, "171210_m37079_mouse2/171210_m37079_mouse2_20171014_790_015na_1hfsds_z5um_1000msec_17-25-29_planes280to299_volume1.tif")]
roi_pths = [os.path.join(src, "171209_f37104_demons/171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13_plane500to519_volume1_v11.RoiSet.zip"),
           os.path.join(src, "171210_m37079_mouse2/171210_m37079_mouse2_20171014_790_015na_1hfsds_z5um_1000msec_17-25-29_planes280to299_volume1_v9.RoiSet.zip")]

#make dict of cell coordinates (IN FIJI TERMS) and profile(s)
pnts_dct = {}


for i, impth in enumerate(impths):
    
    vol = tifffile.imread(impth)
    
    #note python 3 change
    zyx_rois = np.asarray([[int(yy) for yy in xx[0].replace(".roi", "").split("-")] for xx in
                            read_roi_zip(roi_pths[i], include_roi_name=True)])
    
    pnts = np.asarray([[xx[0]-1, xx[1], xx[2]] for xx in zyx_rois])
    
    #save all features of all points in volume
    big_dct = {}
    for pntn in range(len(pnts)):
        
        #save features to dict
        dct = {}
        
        print(pntn)
        try:
            seg = vol[pnts[pntn, 0], pnts[pntn,1]-w:pnts[pntn,1]+w, pnts[pntn,2]-w:pnts[pntn,2]+w]
            assert seg.shape[0] == seg.shape[1]
            cell = np.zeros_like(seg)
            cell[int(cell.shape[0]/2), int(cell.shape[1]/2)] = 100
            
            kernel = np.ones((1,1), np.uint8) 
            cell_dilation = cv2.dilate(cell, kernel, iterations=1) 
        
            #saving
            fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (6,7), gridspec_kw = {"wspace":0, "hspace":0,
                                     "height_ratios": [2,1]})
            
            intensity = vol[pnts[pntn, 0], pnts[pntn,1], pnts[pntn,2]]
            axes[0].imshow(seg*20, "gist_yarg")
            axes[0].imshow(cell_dilation, cmap, alpha = 0.6)
            axes[0].set_ylabel("Location of real cell")
            axes[0].set_title("Volume #{}\nCell coordinate (zyx): {}\nIntensity = {}".format(i+1, tuple((pnts[pntn][0]+1,
                pnts[pntn][1],pnts[pntn][2])), intensity))
            axes[0].set_xticks([], [])
            axes[0].set_yticks([], [])
            
            yprofile = vol[pnts[pntn, 0], pnts[pntn,1]-w:pnts[pntn,1]+w+1, pnts[pntn,2]]
            xprofile = vol[pnts[pntn, 0], pnts[pntn,1], pnts[pntn,2]-w:pnts[pntn,2]+w+1]
            #set z profile padding
            #center z on 0
            center = pnts[pntn, 0] #cell center
            if center >= 11:
                zprofile = vol[center-10:, pnts[pntn,1], pnts[pntn,2]]  
            elif center < 11:
                zprofile = np.ones(xprofile.shape)*float("nan")
                arr = vol[:center+11, pnts[pntn,1], pnts[pntn,2]]  
                zprofile[21-len(arr):] = arr
            
                
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
            
            #save profiles to dct
            dct["yprofile"] = yprofile
            dct["xprofile"] = xprofile
            dct["zprofile"] = zprofile
            dct["cell_center"] = tuple(pnts[pntn])
            dct["intensity"] = intensity
            dct["img_segment"] = seg
            
            big_dct[tuple(pnts[pntn])] = dct
        
        except:
            print("\n")
            print(pnts[pntn])
            print("\n a cell at edge of volume, disregard for now\n\n")
    
    #save points per volume
    pnts_dct[os.path.basename(os.path.dirname(impth))] = big_dct 
    

#write PDF document contains all points
pdf_pages.close()

#save data to pickle
import pickle
sv = "/jukebox/wang/zahra/kelly_cell_detection_analysis/real_cells.p"

with open(sv, "wb") as fp:
    pickle.dump(pnts_dct, fp, protocol=pickle.HIGHEST_PROTOCOL)

data = pickle.load(open(sv, "rb"), encoding = "latin1")