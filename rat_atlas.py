#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 14:57:23 2018

@author: wanglab
"""

import os as sys
from skimage.external import tifffile
import SimpleITK as sitk
import shutil, numpy as np, cv2 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from tools.registration.register import elastix_command_line_call, transformix_command_line_call

if __name__ == '__main__': 

#importing annotations to lightsheet atlas from WHS atlas

    parameters = ['/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/parameterfolder/Order1_Par0000affine.txt',
                  '/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/parameterfolder/Order2_Par0000bspline.txt']
    
    src = '/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/rat_output/median_image_compressed_5.tif'
    atl = '/jukebox/LightSheetData/rat-brody/atlas/modified/WHS_SD_rat_T2star_v1.01_anterior_up_skullremoved_sagittal_CROPPED.tif'
    
    out = '/jukebox/LightSheetData/rat-brody/atlas/meta-atlas/5_brains/mapping'
    elastix_command_line_call(fx = src, mv = atl, out = out, parameters=parameters)



def crop_atlas():

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

