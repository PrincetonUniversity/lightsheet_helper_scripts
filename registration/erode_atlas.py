#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:28:27 2019

@author: wanglab
"""

# Python program to demonstrate erosion and  
# dilation of images. 
import cv2 
import numpy as np 
import SimpleITK as sitk, tifffile
  
# Reading the input image 
img_pth = '/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_zflipped.nrrd'
#img_pth = "/home/wanglab/mounts/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
img = sitk.GetArrayFromImage(sitk.ReadImage(img_pth))
  
# Taking a matrix of size 5 as the kernel 
kernel = np.ones((3,3), np.uint8) 
  
# The first parameter is the original image, 
# kernel is the matrix with which image is  
# convolved and third parameter is the number  
# of iterations, which will determine how much  
# you want to erode/dilate a given image.  
img_erosion = cv2.erode(img, kernel, iterations=3) 

#save out
tifffile.imsave("/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_zflipped_75um_erosion.tif", img_erosion)  
tifffile.imsave("/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_zflipped_75um_erosion.tif", img_erosion)  
