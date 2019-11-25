#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 13:10:18 2019

@author: wanglab
"""

import os, tifffile, cv2, numpy as np, multiprocessing as mp, sys, shutil
from scipy.ndimage import zoom
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.utils.io import load_kwargs
from tools.imageprocessing.orientation import fix_orientation

def downsize(pln, dst):
    
    print(os.path.basename(pln))
    img = tifffile.imread(pln)
    dims = (int(np.ceil(2160/(25/6.5/1.1))), int(np.ceil(4018/(25/6.5/1.1)))) #dims of most of the ones i imaged
    img_resz = cv2.resize(img, dims)
    tifffile.imsave(os.path.join(dst, os.path.basename(pln)), img_resz)
    
def run_downsizing(pth):
    
    print(os.path.basename(pth))
    kwargs = load_kwargs(pth)
    regvol = [xx for xx in kwargs["volumes"] if xx.ch_type == "regch"][0]
    fszdt = regvol.full_sizedatafld_vol
    dst = os.path.join(pth, "downsized_for_atlas") #set destination for downsized planes
    if not os.path.exists(dst): os.mkdir(dst) #make dest directory
    plns = [os.path.join(fszdt, xx) for xx in os.listdir(fszdt) if "tif" in xx]; plns.sort()
    iterlst = [(pln, dst) for pln in plns]
    p = mp.Pool(12)
    p.starmap(downsize, iterlst)

def export_to_tiff(pth, dst):
    fld = os.path.join(pth, "downsized_for_atlas")
    print(os.path.basename(pth))
    plns = [os.path.join(fld, xx) for xx in os.listdir(fld)]; plns.sort()
    arr = np.zeros((len(plns), tifffile.imread(plns[0]).shape[0], tifffile.imread(plns[0]).shape[1]))
    for i,pln in enumerate(plns):
        arr[i] = tifffile.imread(pln)        
        if i%100 == 0: print(i)
    tifffile.imsave(os.path.join(pth, "downsized_for_atlas.tif"), arr.astype("uint16"))
    #remove folder with indiviudal tifs
    shutil.rmtree(fld)
    
    #make 25um iso volume
    arr_dwnsz = zoom(arr, ((2/5), 1, 1), order = 1) #horizontal image, downsizing z by 10um (z step) / 25 um (desired resolution)
    arr_dwnsz_sag = fix_orientation(arr_dwnsz, axes = ("2", "1", "0"))
    tifffile.imsave(os.path.join(dst, os.path.basename(pth)+".tif"), arr_dwnsz_sag.astype("uint16"))
    
    return
    
if __name__ == "__main__":
    
    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    src = "/jukebox//LightSheetData/brodyatlas/processed"
    
    brains = ["k304",
             "c514",
             "w118",
             "c223",
             "k303",
             "k301",
             "e106",
             "k295",
             "w122",
             "h208",
             "k293",
             "k305",
             "k281",
             "k302",
             "k307",
             "c516",
             "h170",
             "k292",
             "e109",
             "c515"]
    
    pths = [os.path.join(src, xx) for xx in brains]

    pth = pths[jobid]
    print(pth)
    #run
    run_downsizing(pth)
    #FIXME: can combine these 2 functions
    dst = "/jukebox/LightSheetData/brodyatlas/atlas/2019_meta_atlas/volumes"
    export_to_tiff(pth, dst)
