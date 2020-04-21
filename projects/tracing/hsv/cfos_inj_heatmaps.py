#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 17:29:23 2020

@author: wanglab
"""

import os, subprocess as sp, tifffile, numpy as np, shutil, matplotlib.pyplot as plt, matplotlib as mpl
from tools.imageprocessing.orientation import fix_orientation
from collections import Counter


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

def find_site(im, thresh=10, filter_kernel=(5,5,5), num_sites_to_keep=1):
    """Find a connected area of high intensity, using a basic filter + threshold + connected components approach
    
    by: bdeverett

    Parameters
    ----------
    img : np.ndarray
        3D stack in which to find site (technically need not be 3D, so long as filter parameter is adjusted accordingly)
    thresh: float
        threshold for site-of-interest intensity, in number of standard deviations above the mean
    filter_kernel: tuple
        kernel for filtering of image before thresholding
    num_sites_to_keep: int, number of injection sites to keep, useful if multiple distinct sites
    
    Returns
    --------
    bool array of volume where coordinates where detected
    """
    from scipy.ndimage.filters import gaussian_filter as gfilt
    from scipy.ndimage import label
    if type(im) == str: im = tifffile.imread(im)

    filtered = gfilt(im, filter_kernel)
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)

    if nlab == 0:
        raise Exception("Site not detected, try a lower threshold?")
    elif nlab == 1:
        return labelled.astype(bool)
    elif num_sites_to_keep == 1:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        return labelled == np.argmax(sizes)+1
    else:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
        return np.in1d(labelled, vals).reshape(labelled.shape)
    
if __name__ == "__main__":
    
    #from tom's jupyter notebook
    thresh=3
    filter_kernel=(3,3,3)
    
    src = "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/injection_site/"
    
    #from tom's code and manual curation... these are the archT+ ones
    
    brains = ["201701_mk01", "201701_mk02", "201701_mk03", "201701_mk04", "201701_mk05",
              "201701_mk06", "201701_mk07", "201701_mk08", "201701_mk10", "201701_mk11"]
    
    #run
    inputlist = [os.path.join(src, xx) for xx in brains]
    
    dct = {"inputlist": inputlist,
      "channel": "01",
      "channel_type": "injch",
      "filter_kernel": (3,3,3), 
      "threshold": 3, 
      "num_sites_to_keep": 1,
      "injectionscale": 45000, 
      "imagescale": 2,
      "reorientation": ("2","0","1"),
      "crop": "[:,:,:]",
      "dst": "/jukebox/wang/zahra/tracing_projects/mapping_paper/figure_data/cfos_heatmap_inj",
      "save_individual": False, 
      "save_tif": True,
      "colormap": "plasma", 
      "atlas": "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif",
      "annotation": "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif",
      "id_table": "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
    }
        
    #segment injection site
    # for brain in inputlist:
    #     vol = os.path.join(os.path.join(src, brain), "transformix_on_inj/result.tif")
    #     arr = find_site(vol, dct["threshold"], dct["filter_kernel"])
    #     tifffile.imsave(os.path.join(dct["dst"], os.path.basename(brain)+".tif"), 
    #                     arr.astype("uint16")*65535) #save to inj site destination
    
    
    #inspect injection sites for the brains i currently have
    imgs = [os.path.join(dct["dst"], xx+".tif") for xx in brains]; imgs.sort()
    #the y-axis cutoff for visualization
    sites = np.array([fix_orientation(tifffile.imread(xx)[:, 420:, :], dct["reorientation"]) for xx in imgs]) 
    atl = fix_orientation(tifffile.imread(dct["atlas"])[:, 420:, :], dct["reorientation"]) #cut at 423 for HSV??
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
        
    limeg = mpl.colors.LinearSegmentedColormap.from_list("", ["white", "limegreen"]) #lime color
    my_cmap = eval("{}(np.arange(plt.cm.Greens.N))".format("limeg"))
    my_cmap[:1,:4] = 0.0  
    my_cmap = mpl.colors.ListedColormap(my_cmap)
    my_cmap.set_under("w")
    plt.figure()
    plt.imshow(np.max(atl, axis=0), cmap="gray")
    plt.imshow(np.max(array, axis=0), alpha=0.90, cmap=my_cmap); plt.colorbar(); plt.axis("off")
    plt.tick_params(length=6)

    plt.savefig("/home/wanglab/Desktop/cfos_inj_heatmap_cb_green.pdf", dpi = 300, transparent = True);
    plt.close()
