#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:42:40 2018

@author: wanglab
"""

#"left" is an approximate

import os

rem_pth = "/jukebox/wang/pisano/tracing_output/antero"
imgd_pth = "/jukebox/wang/pisano/tracing_output/antero_4x"
imgd_prv = "/jukebox/wang/pisano/tracing_output/retro_4x"

img_b4 = [xx for xx in os.listdir(rem_pth)]
img_4x = [yy for yy in os.listdir(imgd_pth)]
img_prv = [ww for ww in os.listdir(imgd_prv)]

left = [zz for zz in img_b4 if zz not in img_4x]

print(len(img_4x))
print(len(img_prv))
print(len(img_b4))
print(len(left))