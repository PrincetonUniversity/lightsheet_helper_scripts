#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 14:57:23 2018

@author: wanglab
"""

import os, numpy as np, sys
from skimage.external import tifffile
from tools.utils.io import *
from tools.registration.register import elastix_command_line_call

if __name__ == '__main__': #do we need this for the commandline?

    parameters = ['/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/parameterfolder/Order1_Par0000affine.txt',
                  '/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/parameterfolder/Order2_Par0000bspline.txt']
    
    src = '/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/rat_output/median_image_compressed_5.tif'
    atl = '/jukebox/LightSheetData/rat-brody/atlas/modified/WHS_SD_rat_T2star_v1.01_anterior_up_skullremoved_sagittal_CROPPED.tif'
    
    out = '/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/mapping'
    elastix_command_line_call(fx = src, mv = atl, out = out, parameters=parameters)