#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:45:02 2019

@author: wanglab
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np, os, sys, subprocess as sp, SimpleITK as sitk, cv2
from skimage.exposure import equalize_adapthist
from scipy.ndimage.filters import gaussian_filter as gfilt
from scipy.ndimage import label
from scipy.ndimage import zoom
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.utils.io import makedir, load_kwargs, listdirfull
from tools.imageprocessing.orientation import fix_orientation

def run_transformix(invol, outpth, transformfile):
    
    #run transformix        
    sp.call(["transformix", "-in", invol, "-out", outpth, "-tp", transformfile])
    print(svlc,"\n   Transformix File Generated: {}".format(outpth))
    
    return outpth

if __name__ == "__main__":
    
    print(sys.argv)
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    src = "/jukebox/wang/pisano/tracing_output/eaat4"
    dst = "/jukebox/wang/zahra/eaat4_screening/201910_analysis"
    
    brains = listdirfull(src)
    
    brain = brains[jobid]
    
    kwargs = load_kwargs(brain)
    cellvol = [vol for vol in kwargs["volumes"] if vol.ch_type == "cellch" or vol.ch_type == "injch"][0]
    fullszfld = cellvol.full_sizedatafld_vol
    
    imgs = [os.path.join(fullszfld, xx) for xx in os.listdir(fullszfld)]; imgs.sort()
    stk = np.array([tif.imread(img) for img in imgs])[:, 1700:, :]
    #stk = tif.imread(src).astype("uint16")
    
    clh = np.array([equalize_adapthist(img, clip_limit=0.05, kernel_size = (50,100), nbins = 65535) for img in stk]) 
    
    power = 3 #how many times to multiple image by itself
    #
    #plt.figure()
    #plt.imshow(clh[300]**power)
    
    thresh = 2
    filter_kernel=(5,5,5)
    
    filtered = gfilt(clh**power, filter_kernel)
    #tif.imsave("/home/wanglab/Desktop/filtered_z200-300.tif", filtered[200:300].astype("float32"))
    #
    #plt.figure()
    #plt.imshow(filtered[300], "gist_yarg")
    
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)
    
    #plt.figure()
    sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
    vals = [i+1 for i in np.argsort(sizes)[::-1]]
    arr = np.in1d(labelled, vals).reshape(labelled.shape)
    
    seg_dst = os.path.join(dst, "rawdata_segmentations"); makedir(seg_dst)
    tif.imsave(os.path.join(seg_dst, "%s.tif" % os.path.basename(brain)), np.in1d(labelled, vals).reshape(labelled.shape).astype("uint16"))
    
    im = tif.imread(imgs[0])
    arr_fullsz = np.zeros((len(imgs), im.shape[0], im.shape[1]))
    print(arr_fullsz.shape)
    print(arr.shape)
    
#################################################################################################################################################################
    arr_fullsz[:, 1700:, :] = arr
    
    #plt.imshow(arr_fullsz[300])
    
    zf,yf,xf = 1,(1/3),(1/3)
    dim = (int(im.shape[1]*xf), int(im.shape[0]*yf))
    arr_resz = np.array([cv2.resize(img, dim) for img in arr_fullsz])
    print(arr_resz.shape)
    
    #to sagittal
    arr_resz_sag = fix_orientation(arr_resz, axes = ("2", "1", "0"))
    
    #now downsample again
    resmpl = tif.imread(cellvol.resampled_for_elastix_vol)
    
    zrsm, yrsm, xrsm = resmpl.shape
    zresz, yresz, xresz = arr_resz_sag.shape
    
    arr_resmpl_sag = zoom(arr_resz_sag, (zrsm/zresz, yrsm/yresz, xrsm/xresz), order = 1)
    
    resz_dst = os.path.join(dst, "downsized_segmentations"); makedir(resz_dst)
    
    tif.imsave(os.path.join(resz_dst, "%s.tif" % os.path.basename(brain)), arr_resmpl_sag.astype("uint16")) 
    
    svlc = os.path.join(brain, "elastix")
    #find transform file
    sig2reg_fld = os.path.join(svlc, os.path.basename(cellvol.downsized_vol))+"/sig_to_reg"
    transformfile = os.path.join(sig2reg_fld, "regtoatlas_TransformParameters.0.txt")
    
    #change the output image type bc otherwise it iterpolates too much and looks weird
    with open(transformfile, "r") as file:
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
    # Write the file out again
    with open(transformfile, "w") as file:
      file.write(filedata)
    
    invol = os.path.join(resz_dst, "%s.tif" % os.path.basename(brain))
    outpth = os.path.join(dst, "transformed_volumes"); makedir(outpth)
    #os.mkdir(outpth)
    outpth = run_transformix(invol, outpth, transformfile)
    
    #cleanup image
    im = sitk.GetArrayFromImage(sitk.ReadImage((os.path.join(outpth, "result.tif"))))
    im[im < 0.1] = 0
    im[im > 0.1] = 255
    tif.imsave(os.path.join(outpth, "%s_trnsfm2atl.tif" % os.path.basename(brain)), im.astype("uint8"))
