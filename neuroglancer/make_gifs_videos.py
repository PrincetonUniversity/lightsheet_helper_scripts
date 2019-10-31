#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:56:34 2019

@author: wanglab
"""

import os, cv2

pth = '/home/wanglab/Documents/neuroglancer/screenshots/20161205_tp_bl6_lob45_1000r_01/amygdala'
list_of_tif_files = [os.path.join(pth, xx) for xx in os.listdir(pth) if "png" in xx]; list_of_tif_files.sort()

#make pngs into video
dst = os.path.join(os.path.dirname(pth), '20161205_tp_bl6_lob45_1000r_01_amygdala.avi')

frame_array = []

fps = 10

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