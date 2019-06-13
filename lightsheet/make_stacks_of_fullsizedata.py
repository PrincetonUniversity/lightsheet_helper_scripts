#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:50:53 2019

@author: wanglab
"""

import SimpleITK as sitk
import tifffile, numpy as np, os

if __name__ == '__main__':
        
    #file that has raw data stacks
    srcs = ["/jukebox/LightSheetData/rat-brody/190523_test_BABB4_2Hr5DCM_4x_488_008na_1hfLs_z5um_100msec_13-56-59",
            "/jukebox/LightSheetData/rat-brody/190523_test_BABB15_1HrDCM_4x_488_008na_1hfLs_z5um_100msec_14-06-55",
            "/jukebox/LightSheetData/rat-brody/190523_test_BABB15_2Hr5DCM_4x_488_008na_1hfLs_z5um_100msec_14-15-46"]
            
    dsts = ["/jukebox/LightSheetData/rat-brody/processed/201906_udisco_tests/BABB4_2Hr5DCM/",
            "/jukebox/LightSheetData/rat-brody/processed/201906_udisco_tests/BABB15_1HrDCM/",
            "/jukebox/LightSheetData/rat-brody/processed/201906_udisco_tests/BABB15_2Hr5DCM/"]

    for src, dst in zip(srcs, dsts):
        
        plns = [os.path.join(src, xx) for xx in os.listdir(src)]; plns.sort()
        
        #get shape    
        y, x = tifffile.imread(plns[3]).shape
        z = len(plns)
        
        if not os.path.exists(dst): os.mkdir(dst)
        
        arr = np.lib.format.open_memmap(os.path.join(dst, "fullsizedata_stack.npy"), dtype = "uint16", mode = "w+", shape = (z,y,x))
            
        #export to memmap array
        for i, pln in zip(range(len(plns)), plns):
            #fill z plns
            print(i)
            img = sitk.GetArrayFromImage(sitk.ReadImage(str(pln)))
            print(img.shape)
            arr[i,:,:] = img
            if i%20 == 0: arr.flush()
        
        #import memmap array to make into tif
        arr = np.lib.format.open_memmap(os.path.join(dst, "fullsizedata_stack.npy"), dtype = "uint16", mode = "r")
        
        #compress and save
        tifffile.imsave(os.path.join(dst, "fullsizedata_stack.tif"), arr, compress = 6)
