#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 13:10:18 2019

@author: wanglab
"""

import os, tifffile, cv2, numpy as np, multiprocessing as mp, sys, shutil
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.utils.io import load_kwargs

src = "/jukebox//LightSheetData/brodyatlas/processed"

brains = ["w118",
         "w122",
         "k293",
         "k302",
         "k307"]

pths = [os.path.join(src, xx) for xx in brains]

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
    plns = [os.path.join(fszdt, xx) for xx in os.listdir(fszdt)]; plns.sort()
    iterlst = [(pln, dst) for pln in plns]
    p = mp.Pool(12)
    p.starmap(downsize, iterlst)

def export_to_tiff(pth):
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
    
if __name__ == "__main__":
    
    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    pth = pths[jobid]
    print(pth)
    #run
    run_downsizing(pth)
    #FIXME: can combine these 2 functions
    export_to_tiff(pth)
    