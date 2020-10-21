#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 09:47:38 2020

@author: wanglab
"""

import tifffile,os
from scipy.ndimage import zoom

#make all volumes the same size?
dst = "/jukebox/LightSheetData/kocher-bee/volume_analysis/volumes_downsized_to_template"

templpth = "/jukebox/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_rotate_croppedZ.tif"
templ = tifffile.imread(templpth)
zt,yt,xt = templ.shape

src = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/experimental_brains"
brs = [os.path.join(src,xx) for xx in os.listdir(src)] #for all brains in the folder

for br in brs:
    print(br)
    img = tifffile.imread(br)
    z,y,x = img.shape
    imgsz = zoom(img, zoom=((zt/z),(yt/y),(xt/x)), order=1)
    tifffile.imsave(os.path.join(dst, os.path.basename(br)), imgsz)

# #read original template
# pth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_rotate.tif"
# img = tifffile.imread(pth)
# #crop in Z
# img_mod = img[:190]
# dst = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_rotate_croppedZ.tif"
# tifffile.imsave(dst, img_mod)

# #do same for annotation
# pth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment.tif"
# img = tifffile.imread(pth)
# #crop in Z
# img_mod = img[:190]
# dst = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment_croppedZ.tif"
# tifffile.imsave(dst, img_mod)

# #read in jacobian
# jcpth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/brain_to_template/Grp16_2.575_elastix/spatialJacobian.tif"
# jac = tifffile.imread(jcpth)

# #%%
# #read in volumes and downsize
# ex1pth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/Grp16_2.575.tif"
# ex2pth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575.tif"
# temppth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_rotate.tif"

# templ = tifffile.imread(temppth)
# templ_downsized = zoom(templ, zoom=(1/3), order=1)
# z,y,x = templ_downsized.shape

# ex1 = tifffile.imread(ex1pth)
# ex1z,ex1y,ex1x = ex1.shape
# ex1_downsized = zoom(ex1, zoom=((z*1.2/ex1z),(y*1.2/ex1y),(x*1.2/ex1x)), order=1)
# ex2 = tifffile.imread(ex2pth)
# ex2z,ex2y,ex2x = ex2.shape
# ex2_downsized = zoom(ex1, zoom=((z*1.2/ex2z),(y*1.2/ex2y),(x*1.2/ex2x)), order=1)

# tifffile.imsave("/home/wanglab/LightSheetData/kocher-bee/volume_analysis/downsized_template_registration/Grp16_2.575.tif", ex1_downsized)
# tifffile.imsave("/home/wanglab/LightSheetData/kocher-bee/volume_analysis/downsized_template_registration/IsoYellow_2.575.tif", ex2_downsized)
# tifffile.imsave("/home/wanglab/LightSheetData/kocher-bee/volume_analysis/downsized_template_registration/Bombus45_2.575umstep_rotate.tif", templ_downsized)

#%%
#moved to reg script
# import sys
# sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
# from tools.registration.register import transformix_command_line_call

# #tranform annotation volume to experimental brain space
# transform = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template_to_brain/Bombus45_2.575umstep_rotate_croppedZ_elastix/TransformParameters.1.txt"

# with open(transform, "r") as file:
#     filedata = file.read()
#     # Replace the target string
#     filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
#     filedata = filedata.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)')
#     # Write the file out again
#     with open(transform, "w") as file:
#       file.write(filedata)
          
# ann = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment_croppedZ.tif"
# dst = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template_to_brain/Bombus45_2.575umstep_rotate_croppedZ_elastix"
# transformix_command_line_call(ann, dst, transform)

#read transformed file
# trsfmpth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template_to_brain/Bombus45_2.575umstep_rotate_croppedZ_elastix/result.tif"
# trsfm = tifffile.imread(trsfmpth)
# #get voxels between 1 and 2 (ideally the curly brain region)
# arr = (trsfm == 2).astype(int)

# #plot
# %matplotlib inline
# import matplotlib.pyplot as plt
# plt.imshow(np.max(arr,axis=0))

# #count total voxels in transform annotation file
# total_voxels_region2_Grp16 = np.sum(arr)
# #compare to original annotation
# total_voxels_region2 = np.sum((tifffile.imread(ann)==2).astype(int))
