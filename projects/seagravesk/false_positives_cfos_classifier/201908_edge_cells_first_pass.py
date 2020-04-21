#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:13:29 2019

@author: wanglab
"""

import scipy.io as sio, os, numpy as np, tifffile, matplotlib.pyplot as plt, cv2, SimpleITK as sitk
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop

"""
Cells look like they are ~6 pixels across (the large ones, anyway), so maybe look at 9 pixels across, with the 
cell center as the center pixel?

Also, maybe makes sense to order the plots by the cell center intensity value, as an approximation of how certain 
we are that what was labeled is actually a cell.

Total cells: 114
"""

#setup
imgs_pth = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/"
mat_pths = [os.path.join(imgs_pth, "cfos/171206_f37077_observer_20171011_790_015na_1hfsds_z5um_1000msec_12-27-06/examples_of_edge_cells_975percentile.mat"),
            os.path.join(imgs_pth, "cfos_201810/181011_f37077_observer_20171011_790_017na_1hfds_z5um_1000msec_13-29-49/examples_of_edge_cells_975percentile.mat"),
            os.path.join(imgs_pth, "cfos/171217_m37109_mouse2_20171018_790_015na_1hfsds_z5um_1000msec_11-51-23/examples_of_edge_cells_975percentile.mat"),
            os.path.join(imgs_pth, "cfos_201810/181023_m37081_observer_20171014_790_017na_1hfds_z5um_1000msec_13-46-16/examples_of_edge_cells_975percentile.mat")
            ]

brain_id = ["DV_f37077_observer_20171011", "VD_f37077_observer_20171011", "DV_m37109_mouse2_20171018", "VD_m37081_observer_20171014"]            

dst = "/jukebox/wang/zahra/kelly_cell_detection_analysis"
pdf_pages = PdfPages(os.path.join(dst, "edge_cell_profiles.pdf")) #compiles into multiple pdfs

#save all volumes

pnts_dct = {}
for i,mat in enumerate(mat_pths):
    dct = sio.loadmat(mat)
    
    arr = dct["edges"]
    
    #find and make volume
    brain = os.path.dirname(mat_pths[i])
    imgs = [os.path.join(brain, xx) for xx in os.listdir(brain) if ".tif" in xx]; imgs.sort()
    
    #only load the 2 planes for which edge cells where selected
    zplns2load = np.unique([xx[0] for xx in arr])
        
    #get all points corresponding to that plane    
    pnts = np.asarray([xx for xx in arr if xx[0] == 700])
    
    #remember layout is y, x (python convention)
    #make dict of cell coordinates (IN FIJI TERMS) and profile(s)
    #save all features of all points in volume
    big_dct = {}
    
    w = 10#window for profiles
    d = 20#depth window
    
    for zpln in zplns2load:
        print(zpln)
        #load plane where point is
        img = tifffile.imread([xx for xx in imgs if "Z{}".format(str(zpln).zfill(4)) in xx])
        
        #make a volume of planes around the point
        depth = 20
        vol = np.zeros((depth*2, img.shape[0], img.shape[1]))
        vol_sz = np.arange(zpln-depth, zpln+depth)
        
        print("\n\n reading image...")
        for z,j in enumerate(vol_sz):
            if z%10 == 0: print(z)
            im = [xx for xx in imgs if "Z{}".format(str(j).zfill(4)) in xx]    
            vol[z] = sitk.GetArrayFromImage(sitk.ReadImage(im))
            
        #get all points corresponding to that plane    
        pnts = np.asarray([[xx[0],xx[1]-1,xx[2]-1] for xx in arr if xx[0] == zpln]) #corrected for matlab 1-based array indices
    
        for pntn in pnts:
            
            print(pntn)
            
            #save features to dict
            dct = {}
            z,y,x = pntn
            seg = img[y-d:y+d, x-d:x+d]
            cell = np.zeros_like(seg)
            cell[int(cell.shape[0]/2), int(cell.shape[1]/2)] = 100
            
            kernel = np.ones((1,1), np.uint8) 
            cell_dilation = cv2.dilate(cell, kernel, iterations=1) 
        
            #saving
            fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (6,7), gridspec_kw = {"wspace":0, "hspace":0,
                                     "height_ratios": [2,1]})
            
            axes[0].imshow(seg*20, "gist_yarg")
            axes[0].imshow(cell_dilation, cmap, alpha = 0.6)
            axes[0].set_ylabel("Location of edge cell")
            axes[0].set_title("Cell coordinate (zyx): {}".format(tuple(pntn)))
            axes[0].set_xticks([], [])
            axes[0].set_yticks([], [])
            
            intensity = img[y, z]
            yprofile = img[y-w:y+w+1, x]
            xprofile = img[y, x-w:x+w+1]
            zprofile = vol[d-w:d+w+1, y, x]
            
            axes[1].plot(zprofile, label = "Z")
            axes[1].plot(yprofile, label = "Y")
            axes[1].plot(xprofile, label = "X")
            axes[1].set_ylabel("Pixel intensity")
            axes[1].set_xlabel("Distance (pixels)")
            axes[1].set_xticks(np.arange(0, len(xprofile), 2))
            axes[1].legend()
            
            #done with the page
            pdf_pages.savefig(dpi = 300, bbox_inches = 'tight') 
            plt.close()
            
            #save profiles to dct
            dct["yprofile"] = yprofile
            dct["xprofile"] = xprofile
            dct["zprofile"] = zprofile
            dct["intensity"] = intensity
            dct["cell_center"] = tuple(pntn)
            dct["img_segment"] = seg
            
            big_dct[tuple(pntn)] = dct
            
    pnts_dct[brain_id[i]] = big_dct 
    
#write PDF document contains all points
pdf_pages.close()

print("\nexporting to pickle...\n")
#save data to pickle
import pickle
sv = os.path.join(dst, "edge_cells.p")

with open(sv, "wb") as fp:
    pickle.dump(pnts_dct, fp, protocol=pickle.HIGHEST_PROTOCOL)

#data = pickle.load(open(sv, "rb"), encoding = "latin1")
