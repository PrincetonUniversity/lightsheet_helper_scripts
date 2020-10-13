#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 17:37:06 2020

@author: wanglab
"""

import os, subprocess as sp, tifffile, numpy as np, shutil, matplotlib.pyplot as plt, matplotlib as mpl
from tools.utils.io import makedir, load_kwargs, listdirfull
from tools.imageprocessing.orientation import fix_orientation
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

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
    
def run_transformix(invol, outpth, transformfile):
    
    #run transformix        
    sp.call(["transformix", "-in", invol, "-out", outpth, "-tp", transformfile])
    print("\n   Transformix File Generated: {}".format(outpth))
    
    return outpth

#here, im using the injection detection code for segmenting out the injection site from the resampled volume - more accurate
#later i can use that segmented volume in the downsampled space and convert it to atlas space

if __name__ == "__main__":


    #run
    inputlist = [
        "/jukebox/LightSheetTransfer/tp/20200930_17_32_58_hsv_36hr_7/Ex_642_Em_2/downsized_for_atlas.tif",
        "/jukebox/LightSheetTransfer/tp/20201001_10_57_49_hsv_36h_6/Ex_642_Em_2/downsized_for_atlas.tif",
        "/jukebox/LightSheetTransfer/tp/20201001_17_13_35_hsv_28h_2/Ex_642_Em_2/downsized_for_atlas.tif",
        "/jukebox/LightSheetTransfer/tp/20201001_15_39_26_hsv_28h_4/Ex_642_Em_2/downsized_for_atlas.tif",
        "/jukebox/LightSheetTransfer/tp/PRV_50hr-019/Ex_642_Em_2/downsized_for_atlas.tif",
        "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_PRV/natneuroreviews_tompisano_PRV_36hr-015/imaging_request_1/output/processing_request_1/resolution_4x/Ex_642_Em_2/downsized_for_atlas.tif",
        "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_PRV/natneuroreviews_tompisano_PRV_28hr-011/imaging_request_1/output/processing_request_1/resolution_4x/Ex_642_Em_2/downsized_for_atlas.tif",
        "/jukebox/wang/pisano/tracing_output/bl6_ts/20150804_tp_bl6_ts04/20150804_tp_bl6_ts04_555_z3um_70msec_3hfds_resized_ch00_resampledforelastix.tif",
        "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_CTB/natneuroreviews_tompisano_CTB-001/imaging_request_1/output/processing_request_1/resolution_4x/Ex_561_Em_1/downsized_for_atlas.tif",
        "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_CTB/natneuroreviews_tompisano_CTB-002/imaging_request_1/output/processing_request_1/resolution_4x/Ex_561_Em_1/downsized_for_atlas.tif"
        ]
    
    transformfiles = [
        "/jukebox/LightSheetTransfer/tp/20200930_17_32_58_hsv_36hr_7/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/LightSheetTransfer/tp/20201001_10_57_49_hsv_36h_6/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/LightSheetTransfer/tp/20201001_17_13_35_hsv_28h_2/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/LightSheetTransfer/tp/20201001_15_39_26_hsv_28h_4/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/LightSheetTransfer/tp/PRV_50hr-019/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_PRV/natneuroreviews_tompisano_PRV_36hr-015/imaging_request_1/output/processing_request_1/resolution_4x/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_PRV/natneuroreviews_tompisano_PRV_28hr-011/imaging_request_1/output/processing_request_1/resolution_4x/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/wang/pisano/tracing_output/bl6_ts/20150804_tp_bl6_ts04/elastix_inverse_transform/injch_20150804_tp_bl6_ts04_488w_647_z3um_250msec_1hfds/20150804_tp_bl6_ts04_488w_647_z3um_250msec_1hfds_resized_ch01_resampledforelastix_atlas2reg2sig/reg2sig_TransformParameters.1.txt",
        "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_CTB/natneuroreviews_tompisano_CTB-001/imaging_request_1/output/processing_request_1/resolution_4x/elastix_inverse_transform/TransformParameters.1.txt",
         "/jukebox/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_CTB/natneuroreviews_tompisano_CTB-002/imaging_request_1/output/processing_request_1/resolution_4x/elastix_inverse_transform/TransformParameters.1.txt"
         ]
        
    brains = ["hsv_36hr_7","hsv_36h_6",
              "hsv_28h_2","hsv_28h_4","PRV_50hr-019","PRV_36hr-015",
              "PRV_28hr-011", "20150804_tp_bl6_ts04", "CTB-001", "CTB-002"]
    
    dct = {"inputlist": inputlist,
      "filter_kernel": (5,5,5), 
      "threshold": 10, 
      "injectionscale": 45000, 
      "imagescale": 2,
      "reorientation": ("2","0","1"),
      "crop": "[:,450:,:]",
      "dst": "/home/wanglab/Desktop/short_time_series_inj",
      "save_individual": True, 
      "save_tif": True,
      "colormap": "plasma", 
      "atlas": "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif",
      "annotation": "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif",
      "id_table": "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"
    }
        
    #segment injection site
    for i, vol in enumerate(inputlist):
        arr = find_site(vol, dct["threshold"], dct["filter_kernel"])
        tifffile.imsave(os.path.join(dct["dst"], brains[i]+".tif"), arr.astype("uint16")*65535) #save to inj site destination
    
    #transform atlas to segmentation space
    #this is because i did not run the reg --> atlas transform
    #btw this will take long, and i don't recommend parallelizing bc of transformix
    for i, vol in enumerate(inputlist):
        invol = dct["annotation"]
        #run transformix
        outpth = os.path.join(dct["dst"], brains[i]); makedir(outpth)
        outpth = run_transformix(invol, outpth, transformfiles[i])
        
       #save out
        tifffile.imsave(os.path.join(dct["dst"],brains[i]+"_annotation.tif"),
                        img.astype("uint16"))
        shutil.rmtree(outpth) #delete original transformed file
    
    
#    #inspect injection sites for the brains i currently have
#    imgs = listdirfull(dct["dst"], "tif"); imgs.sort()
#    sites = np.array([fix_orientation(tifffile.imread(xx)[:, 450:, :], dct["reorientation"]) for xx in imgs]) #the y-axis cutoff for visualization
#    atl = fix_orientation(tifffile.imread(dct["atlas"])[:, 450:, :], dct["reorientation"])
#    
#    my_cmap = eval("plt.cm.{}(np.arange(plt.cm.RdBu.N))".format("viridis"))
#    my_cmap[:1,:4] = 0.0  
#    my_cmap = mpl.colors.ListedColormap(my_cmap)
#    my_cmap.set_under("w")
#    plt.figure()
#    plt.imshow(np.max(atl, axis=0), cmap="gray")
#    plt.imshow(np.max(np.sum(sites, axis=0), axis = 0), alpha=0.90, cmap=my_cmap); plt.colorbar(); plt.axis("off")
#    
#    plt.savefig(os.path.join("/home/wanglab/Desktop/heatmap.pdf"), dpi = 300, transparent = True);
#    plt.close()