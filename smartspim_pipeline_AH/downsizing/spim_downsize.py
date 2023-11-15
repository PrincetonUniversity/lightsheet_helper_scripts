#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:04:02 2020

@author: wanglab
"""

import os, numpy as np, tifffile as tif, SimpleITK as sitk, cv2, multiprocessing as mp, shutil, sys
from scipy.ndimage import zoom

def fast_scandir(dirname):
    """ gets all folders recursively """
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders

def resize_helper(img, dst, resizef):
    print(os.path.basename(img))
    im = sitk.GetArrayFromImage(sitk.ReadImage(img))
    y,x = im.shape
    yr = int(y/resizef); xr = int(x/resizef)
    im = cv2.resize(im, (xr, yr), interpolation=cv2.INTER_LINEAR)
    tif.imsave(os.path.join(dst, os.path.basename(img)), 
                    im.astype("uint16"), compress=1)

if __name__ == "__main__":
    
    #takes 1 command line args
    print(sys.argv)
    src=str(sys.argv[1]) #folder to main image folder
    try:
        pth = fast_scandir(src)[-1] #gets to the end of directory tree
    except:
        pth = src #if images already in first directory
    print("\nPath to stitched images: %s\n" % pth)
    #path to store downsized images
    dst = str(sys.argv[2])
    print("\nPath to storage directory: %s\n\n" % dst)
    if not os.path.exists(dst): os.mkdir(dst)
    imgs = [os.path.join(pth, xx) for xx in os.listdir(pth) if "tif" in xx]
    z = len(imgs)
    resizef = 5 #factor to downsize imgs by
    iterlst = [(img, dst, resizef) for img in imgs]
    p = mp.Pool(12)
    p.starmap(resize_helper, iterlst)
    p.terminate()
    
    #now downsample to 140% of pma atlas
    imgs = [os.path.join(dst, xx) for xx in os.listdir(dst) if "tif" in xx] 
    imgs.sort() 
    dv = str(sys.argv[3])
    if (dv == 'v'): 
        imgs.reverse()
    
    z = len(imgs)
    y,x = sitk.GetArrayFromImage(sitk.ReadImage(imgs[0])).shape
    arr = np.zeros((z,y,x))

    atl_name=str(sys.argv[4])
    atlpth = ""
    if (atl_name == "PMA"):
        atlpth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
    elif (atl_name == "Allen"):
        atlpth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
    elif (atl_name == "cb"):
        atlpth = "/jukebox/LightSheetTransfer/atlas/cb_sagittal_atlas_20um_iso.tif"
    elif (atl_name == "PRA"):
        atlpth = "/jukebox/brody/lightsheet/atlasdir/mPRA.tif"
    elif (atl_name == "cz"):
        atlpth = "/jukebox/witten/Chris/data/clearmap2/utilities/allen-atlas-cz/average_template_25_sagittal_forDVscans_cz.tif"
    elif (atl_name == "hem"):
        atlpth = "/jukebox/LightSheetTransfer/hem_sagittal_atlas.tif"
    else:
        raise ValueError("Specified atlas does not exist")

    atl = sitk.GetArrayFromImage(sitk.ReadImage(atlpth))
    atlz,atly,atlx = atl.shape #get shape, sagittal
    #read all the downsized images
    for i,img in enumerate(imgs):
        if i%5000==0: print(i)
        arr[i,:,:] = sitk.GetArrayFromImage(sitk.ReadImage(img)) #horizontal
    #switch to sagittal
    arrsag = np.swapaxes(arr,2,0)
    z,y,x = arrsag.shape
    print((z,y,x))
    print("\n**********downsizing....heavy!**********\n")
    
    arrsagd = zoom(arrsag, ((atlz*1.4/z),(atly*1.4/y),(atlx*1.4/x)), order=1)
    shutil.rmtree(dst)
    os.mkdir(dst)
    tif.imsave(os.path.join(dst, "downsized_for_atlas.tif"), arrsagd.astype("uint16"))

    # print("\ndeleting storage directory after making volume...\n %s" % dst)
