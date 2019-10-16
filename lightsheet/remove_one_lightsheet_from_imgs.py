#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 17:47:32 2019

@author: wanglab
"""

import os, shutil

pth = "/jukebox/LightSheetTransfer/Jess/opto_ai27d"
dst = "/jukebox/LightSheetTransfer/Jess/opto_ai27d_right_sheet"

if not os.path.exists(dst): os.mkdir(dst)

pths = [os.path.join(pth, xx) for xx in os.listdir(pth)]

for brain in pths:
    print(brain)
    fls_to_mv = [os.path.join(brain, xx) for xx in os.listdir(brain) if "C01" in xx]
    
    brain_dst = os.path.join(dst, os.path.basename(brain))
    if not os.path.exists(brain_dst): os.mkdir(brain_dst)
    
    for fl in fls_to_mv:
        shutil.move(fl, brain_dst)
