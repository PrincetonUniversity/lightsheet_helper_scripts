#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:30:20 2020

@author: wanglab
"""

import os, numpy as np, tifffile as tif, SimpleITK as sitk, cv2, multiprocessing as mp, shutil, sys
from scipy.ndimage import zoom
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.registration.register import elastix_command_line_call

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
    
    #array job paralellization across brains
    src = "/jukebox/LightSheetData/wang-mouse/seagravesk"
    #channel to downsize
    ch = "Ex_785_Em_3"
    brains = [os.path.join(src,xx) for xx in os.listdir(src) if os.path.exists(os.path.join(src,xx,ch,"stitched"))]
    brains.sort()
    #paralellize across brains
    # print(os.environ["SLURM_ARRAY_TASK_ID"])
    # jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    #select brain to run
    brain = os.path.join(src, "20200916_19_25_35_f37080_mouse1_20171015")#brains[jobid]
    #use corrected images if available
    try:
        pth = fast_scandir(os.path.join(brain,ch,"corrected"))[-1]
    except:
        pth = fast_scandir(os.path.join(brain,ch,"stitched"))[-1]
    print("\npath to stitched images: %s\n" % pth)
    #path to store downsized images
    dst = os.path.join(brain, ch, "downsized")
    print("\npath to storage directory: %s\n\n" % dst)
    if not os.path.exists(dst): os.mkdir(dst)
    imgs = [os.path.join(pth, xx) for xx in os.listdir(pth) if "tif" in xx]
    z = len(imgs)
    resizef = 5 #factor to downsize imgs by
    iterlst = [(img, dst, resizef) for img in imgs]
    p = mp.Pool(12)
    p.starmap(resize_helper, iterlst)
    p.terminate()
    #now downsample to reg volume imaged previously (reg to allen atlas)
    imgs = [os.path.join(dst, xx) for xx in os.listdir(dst) if "tif" in xx]; imgs.sort()
    z = len(imgs)
    y,x = sitk.GetArrayFromImage(sitk.ReadImage(imgs[0])).shape
    arr = np.zeros((z,y,x))
    #get dimensions of registration volume from previously imaging run; these are the same for all the volumes
    regz,regy,regx = 592,686,416 #get shape, sagittal
    #read all the downsized images
    for i,img in enumerate(imgs):
        if i%5000==0: print(i)
        arr[i,:,:] = sitk.GetArrayFromImage(sitk.ReadImage(img)) #horizontal
    #switch to sagittal
    arrsag = np.swapaxes(arr,2,0)
    z,y,x = arrsag.shape
    print((z,y,x))
    print("\n**********downsizing....heavy!**********\n")
    arrsagd = zoom(arrsag, ((regz/z),(regy/y),(regx/x)), order=1)
    #note, custom path
    tif.imsave(os.path.join(os.path.dirname(dst), "%s_downsized_for_atlas.tif" % ch), arrsagd.astype("uint16"))
    print("\ndeleting storage directory after making volume...\n %s" % dst)
    shutil.rmtree(dst)
    #run registration of 790 channel to old reg channel volume
    print("\n**********running registration...**********\n\n")
    fx = os.path.join(os.path.dirname(dst), "%s_downsized_for_atlas.tif" % ch)
    volsrc = "/jukebox/wang/seagravesk/lightsheet/201710_cfos_left_side_only_registration"
    brnm = os.path.basename(brain)[18:31] #remove date and name tag
    br = os.path.join(volsrc,brnm)
    mv = [os.path.join(br,xx) for xx in os.listdir(br) if "488" in xx and "resampled" in xx][0]
    #copy 488 volume to current folder 
    print("\n**********copying old reg volume to new processed directory**********\n\n")
    shutil.copy(mv, os.path.join(os.path.dirname(dst), "reg_downsized_for_atlas.tif"))
    print("\npath to downsized vol for inverse registration to reg channel: %s" % fx)
    print("\npath to reg channel: %s" % mv)
    out = os.path.join(brain, "elastix_inverse_transform")
    if not os.path.exists(out): os.mkdir(out)
    #set params
    param_fld = "/jukebox/LightSheetData/wang-mouse/seagravesk/20200810_13_10_58_f37080_mouse2_20171015_slow_focus/parameters" #change if using rat
    params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    #run
    e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
