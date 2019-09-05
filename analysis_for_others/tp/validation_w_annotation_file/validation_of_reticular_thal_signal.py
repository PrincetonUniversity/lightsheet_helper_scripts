#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 17:06:00 2019

@author: wanglab
"""

import tifffile as tif, numpy as np

ann = tif.imread("/home/wanglab/mounts/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")

rtn = ann == 262
vpm = ann == 733
vpl = ann == 718

tif.imsave("/home/wanglab/Desktop/rtn.tif", rtn.astype("uint16"))
tif.imsave("/home/wanglab/Desktop/vpl.tif", vpl.astype("uint16"))
tif.imsave("/home/wanglab/Desktop/vpm.tif", vpm.astype("uint16"))

