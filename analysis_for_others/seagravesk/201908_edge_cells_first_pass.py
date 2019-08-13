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
(hand-picked by Kelly on 20190807, based on the automatically-labeled cells in the file called 
‘sliding_diff_peak_find_975percentile_20190227_format2.mat’)
Number of cells: 26
Contains a (z,y,x) matrix called ‘edges’.

Cells look like they are ~6 pixels across (the large ones, anyway), so maybe look at 9 pixels across, with the 
cell center as the center pixel?

Also, maybe makes sense to order the plots by the cell center intensity value, as an approximation of how certain 
we are that what was labeled is actually a cell.

"""

#setup
imgs_pth = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/"
mat = os.path.join(imgs_pth, "171206_f37077_observer_20171011_790_015na_1hfsds_z5um_1000msec_12-27-06/examples_of_edge_cells_975percentile.mat")

dct = sio.loadmat(mat)

arr = dct["edges"]

brainname = "171206_f37077_observer_20171011_790_015na_1hfsds_z5um_1000msec_12-27-06"
brain = os.path.join(imgs_pth, brainname)

imgs = [os.path.join(brain, xx) for xx in os.listdir(brain) if ".tif" in xx]; imgs.sort()

#only load the 2 planes for which edge cells where selected
zplns2load = np.unique([xx[0] for xx in arr])


#zpln 700
zpln = zplns2load[0]
#load plane where point is
img = tifffile.imread([xx for xx in imgs if "Z{}".format(str(zplns2load[0]).zfill(4)) in xx])

#make a volume of planes around the point
depth = 50
vol = np.zeros((depth*2, img.shape[0], img.shape[1]))
vol_sz = np.arange(zpln-depth, zpln+depth)

for z,i in enumerate(vol_sz):
    if z%10 == 0: print(z)
    im = [xx for xx in imgs if "Z{}".format(str(i).zfill(4)) in xx]    
    vol[z] = sitk.GetArrayFromImage(sitk.ReadImage(im))
    
#get all points corresponding to that plane    
pnts = np.asarray([xx for xx in arr if xx[0] == 700])

#determine the window around which you want to look (pixels)
#remember layout is y, x (python convention)
#y = 0, x = 1
#testing
#pntn = 11
#w = 50
#seg = img[pnts[pntn,1]-w:pnts[pntn,1]+w, pnts[pntn,2]-w:pnts[pntn,2]+w]
#cell = np.zeros_like(seg)
#cell[int(cell.shape[0]/2), int(cell.shape[1]/2)] = 100
#
#kernel = np.ones((3,3), np.uint8) 
#cell_dilation = cv2.dilate(cell, kernel, iterations=1) 
#
#plt.imshow(seg*10, "gist_yarg")
#plt.imshow(cell_dilation, cmap, alpha = 0.3)
#plt.axis("off")
#
#plt.figure()
#yprofile = img[pnts[pntn,1]-w:pnts[pntn,1]+w, pnts[pntn,2]]
#xprofile = img[pnts[pntn,1], pnts[pntn,2]-w:pnts[pntn,2]+w]
#zprofile = vol[:, pnts[pntn,1], pnts[pntn,2]]
#plt.plot(yprofile)
#plt.plot(xprofile)
#plt.plot(zprofile)
#plt.close()

#make dict of cell coordinates (IN FIJI TERMS) and profile(s)
pnts_dct = {}
pnts_dct["brainname"] = brainname

w = 10
d = 20
dst = "/home/wanglab/Desktop/edge_cell_profiles.pdf"
pdf_pages = PdfPages(dst) #compiles into multiple pdfs

for zpln in zplns2load:
    print(zpln)
    #load plane where point is
    img = tifffile.imread([xx for xx in imgs if "Z{}".format(str(zpln).zfill(4)) in xx])
    
    #make a volume of planes around the point
    depth = 20
    vol = np.zeros((depth*2, img.shape[0], img.shape[1]))
    vol_sz = np.arange(zpln-depth, zpln+depth)
    
    print("\n\n reading image...")
    for z,i in enumerate(vol_sz):
        if z%10 == 0: print(z)
        im = [xx for xx in imgs if "Z{}".format(str(i).zfill(4)) in xx]    
        vol[z] = sitk.GetArrayFromImage(sitk.ReadImage(im))
        
    #get all points corresponding to that plane    
    pnts = np.asarray([[xx[0],xx[1]-1,xx[2]-1] for xx in arr if xx[0] == zpln]) #corrected for matlab 1-based array indices

    for pntn in range(len(pnts)):
        
        print(pntn)
        
        #save features to dict
        dct = {}
        
        seg = img[pnts[pntn,1]-d:pnts[pntn,1]+d, pnts[pntn,2]-d:pnts[pntn,2]+d]
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
        axes[0].set_title("Cell coordinate (zyx): {}".format(tuple((pnts[pntn][0],pnts[pntn][1],pnts[pntn][2]))))
        axes[0].set_xticks([], [])
        axes[0].set_yticks([], [])
        
        intensity = img[pnts[pntn,1], pnts[pntn,2]]
        yprofile = img[pnts[pntn,1]-w:pnts[pntn,1]+w, pnts[pntn,2]]
        xprofile = img[pnts[pntn,1], pnts[pntn,2]-w:pnts[pntn,2]+w]
        zprofile = vol[d-w:d+w, pnts[pntn,1], pnts[pntn,2]]
        
        axes[1].plot(zprofile, label = "Z")
        axes[1].plot(yprofile, label = "Y")
        axes[1].plot(xprofile, label = "X")
        axes[1].set_ylabel("Pixel intensity")
        axes[1].set_xlabel("Distance (pixels)")
        axes[1].set_xticks(np.arange(0, len(xprofile), 2))
        axes[1].legend()
        
        #done with the page
        pdf_pages.savefig(dpi = 300, bbox_inches = 'tight') 
        
        #save profiles to dct
        dct["yprofile"] = yprofile
        dct["xprofile"] = xprofile
        dct["zprofile"] = zprofile
        dct["intensity"] = intensity
        dct["cell_center"] = tuple(pnts[pntn])
        dct["img_segment"] = seg
        
        pnts_dct[tuple(pnts[pntn])] = dct
        
#write PDF document contains all points
pdf_pages.close()

#save data to pickle
import pickle
sv = "/home/wanglab/Desktop/edge_cells.p"

with open(sv, "wb") as fp:
    pickle.dump(pnts_dct, fp, protocol=pickle.HIGHEST_PROTOCOL)

#data = pickle.load(open(sv, "rb"), encoding = "latin1")
