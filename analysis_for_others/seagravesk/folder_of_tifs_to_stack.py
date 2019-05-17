#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 11:53:51 2019

@author: wanglab
"""

import tifffile, os, numpy as np

fld = "/home/wanglab/mounts/wang/seagravesk/lightsheet/testing_smartSPIM_20190506/m57202_4x/642"

fls = [os.path.join(fld, xx) for xx in os.listdir(fld)]
fls.sort()

shape = (len(fls), tifffile.imread(fls[0]).shape[0], tifffile.imread(fls[0]).shape[1])

#init 
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/m57202_smartSPIM_4x.npy", dtype = "uint16", mode = "w+", shape = shape)

for i, fl in enumerate(fls):
    print(os.path.basename(fl))
    if i > 952:
        slc = tifffile.imread(fl)
        arr[i] = slc
        if i%1 == 0: arr.flush()