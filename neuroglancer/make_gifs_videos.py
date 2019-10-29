#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:56:34 2019

@author: wanglab
"""

import os, cv2

pth = '/home/wanglab/Documents/neuroglancer/screenshots/20170411_db_bl6_crii_mid_53hr/rtn_vpm_vpl_md'
list_of_tif_files = [os.path.join(pth, xx) for xx in os.listdir(pth) if "png" in xx]; list_of_tif_files.sort()

#make pngs into video
dst = os.path.join(pth, '20170411_db_bl6_crii_mid_53hr_rtn_vpm_vpl_md.avi')

frame_array = []

fps = 20

for png in list_of_tif_files:
    #reading each files
    img = cv2.imread(png)
    height, width, layers = img.shape
    size = (width,height)
    
    #inserting the frames into an image array
    frame_array.append(img)

out = cv2.VideoWriter(dst, cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
    
for i in range(len(frame_array)):
    # writing to a image array
    out.write(frame_array[i])
out.release()