#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 12:47:31 2019

@author: wanglab
"""

import numpy as np
import tifffile

arr_pth = "/jukebox/scratch/zmd/m57202_smartSPIM_4x.npy"

arr = np.lib.format.open_memmap(arr_pth, dtype = "float32", mode = "r")

tifffile.imsave("/jukebox/scratch/zmd/m57202_smartSPIM_4x.tif", arr, compress = 6)