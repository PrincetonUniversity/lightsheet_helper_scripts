#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 13:49:34 2019

@author: wanglab
"""

import os, tifffile, cv2, numpy as np, multiprocessing as mp, sys, shutil
from scipy.ndimage import zoom
sys.path.append('/jukebox/wang/zahra/python/lightsheet_py3')
from tools.utils.io import load_kwargs
from tools.imageprocessing.orientation import fix_orientation

src = '/jukebox//LightSheetData/brodyatlas/processed'

brains = ['w118',
         'w122',
         'k293',
         'k302',
         'k307']

inputs = [os.path.join(src, xx+"/downsized_for_atlas.tif") for xx in brains]

output_fld = '/jukebox//LightSheetData/brodyatlas/atlas/2019_meta_atlas'
if not os.path.exists(output_fld): os.mkdir(output_fld)

data_fld = '/jukebox//LightSheetData/brodyatlas/atlas/2019_meta_atlas/volumes'
if not os.path.exists(data_fld): os.mkdir(data_fld)

for vol in inputs:
    arr = tifffile.imread(vol)
    arr_dwnsz = zoom(arr, ((2/5), 1, 1), order = 1) #horizontal image, downsizing z by 10um (z step) / 25 um (desired resolution)
    arr_dwnsz_sag = fix_orientation(arr_dwnsz, axes = ("2", "1", "0"))
    tifffile.imsave(os.path.join(data_fld, os.path.basename(os.path.dirname(vol))+".tif"), arr_dwnsz_sag)
#%%
#registration to seed
parameters = '/jukebox/wang/zahra/python/lightsheet_py3/parameterfolder' #start with basic affine/bspile