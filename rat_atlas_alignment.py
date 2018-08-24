#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 13:30:28 2018

@author: wanglab 
"""

import SimpleITK as sitk
import shutil, numpy as np, cv2 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from tools.registration.register import elastix_command_line_call, transformix_command_line_call

# the way i tried it....
#sitkvol=sitk.ReadImage( '/home/wanglab/mounts/wang/zahra/zahra_lightsheet_copy/WHS_SD_rat_atlas_v2_anterior_up_sagittal_cropped.tif')
#sitk.Show(vol,'atlas',debugOn=True)
#vol = sitk.GetArrayFromImage(sitkvol) <--- makes images into numpy array

# tom's way
#load atlas file
pth = '/home/wanglab/mounts/LightSheetData/rat-brody/atlas/modified/WHS_SD_rat_T2star_v1.01_anterior_up_skullremoved_sagittal.tif'
from skimage.external import tifffile
vol = tifffile.imread(pth) #<--- tifffile converts tif into numpy array
#sitk.Show(vol,'atlas')

vol = vol[:, 26:800, 62:364] #crops the atlas image in the approximate dimensions you would like
#order of dimensions is z,y,x
tifffile.imsave('/home/wanglab/mounts/wang/zahra/example_elastix_script/data/WHS_SD_rat_T2star_v1.01_anterior_up_skullremoved_sagittal_CROPPED_black_border.tif', vol) #saves it in the directory