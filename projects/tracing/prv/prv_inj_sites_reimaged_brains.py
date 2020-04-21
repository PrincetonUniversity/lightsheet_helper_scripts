#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 17:32:54 2020

@author: wanglab
"""

import os, subprocess as sp, tifffile, numpy as np, shutil, matplotlib.pyplot as plt, matplotlib as mpl
from tools.analysis.analyze_injection_inverse_transform import pool_injections_inversetransform
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

    src = "/jukebox/wang/zahra/tracing_projects/prv/prv_inj_vols_from_reimaged_brains"
    
    brains = ["20180215_jg_bl6f_prv_08"]
    
    #run
    inputlist = [os.path.join(src, xx) for xx in brains]
    
    dct = {"inputlist": inputlist,
      "channel": "01",
      "channel_type": "injch",
      "filter_kernel": (3,3,3), 
      "threshold": 8, 
      "num_sites_to_keep": 1,
      "injectionscale": 45000, 
      "imagescale": 2,
      "reorientation": ("2","0","1"),
      "crop": "[:,:,:]",
      "dst": "/jukebox/wang/zahra/tracing_projects/prv/prv_injection_sites",
      "save_individual": False, 
      "save_tif": True,
      "colormap": "plasma", 
      "atlas": "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif",
      "annotation": "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif",
      "id_table": "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"
    }
        
    #segment injection site
    for brain in inputlist:
        vol = [os.path.join(brain, xx) for xx in os.listdir(brain) if "tif" in xx][0]
        arr = find_site(vol, dct["threshold"], dct["filter_kernel"])
        tifffile.imsave(os.path.join(dct["dst"], os.path.basename(brain)+".tif"), arr.astype("uint16")*65535) #save to inj site destination
    
    #even though we have the voxel counts in the csv file, i would still prefer to have the registered volumes just in case
    #that is how we did the segmentation for the h129 anyways
    #btw this will take long, and i don't recommend parallelizing bc of transformix
    for fld in inputlist:
        
        svlc = [os.path.join(fld, xx) for xx in os.listdir(fld) if "tif" not in xx][0]
        #find transform file
        sig2reg_fld = os.path.join(svlc, "sig_to_reg")
        transformfiles = [os.path.join(sig2reg_fld, "TransformParameters.0.txt"),
                          os.path.join(sig2reg_fld, "TransformParameters.1.txt"),
                          os.path.join(sig2reg_fld, "regtoatlas_TransformParameters.0.txt"),
                          os.path.join(sig2reg_fld, "regtoatlas_TransformParameters.1.txt")]
        
        for transformfile in transformfiles:
            #change the output image type bc otherwise it iterpolates too much and looks weird
            with open(transformfile, "r") as file:
                filedata = file.read()
            # Replace the target string
            filedata = filedata.replace("/jukebox/wang/pisano/tracing_output/retro_4x/"+os.path.basename(fld)+"/elastix", fld)
            filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
            # Write the file out again
            with open(transformfile, "w") as file:
              file.write(filedata)
        
        invol = os.path.join(dct["dst"], os.path.basename(fld)+".tif")
        #run inj detection
        outpth = os.path.join(dct["dst"], os.path.basename(fld)); makedir(outpth)
        outpth = run_transformix(invol, outpth, os.path.join(sig2reg_fld, "regtoatlas_TransformParameters.1.txt"))
        
        #fix the negative #'s around the site, and overwrite the tif onto the resampled version
        #this was checked and is consistent w the transform
        img = tifffile.imread(outpth+"/result.tif")
        img[img < 0] = 0
        assert np.sum(img < 0) == 0 #make sure there are no negative values
        tifffile.imsave(invol, img.astype("uint16"))
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