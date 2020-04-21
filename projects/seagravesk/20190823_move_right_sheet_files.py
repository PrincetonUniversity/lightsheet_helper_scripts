#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 14:39:58 2019

@author: wanglab
"""

import os, shutil

pth = "/jukebox/LightSheetTransfer/kelly/201908_cfos"
dst = "/jukebox/LightSheetTransfer/kelly/201908_cfos_right_sheet"
if not os.path.exists(dst): os.mkdir(dst)

pths = [os.path.join(pth, xx) for xx in os.listdir(pth)]

for brain in pths:
    fls_to_mv = [os.path.join(brain, xx) for xx in os.listdir(brain) if "C01" in xx]
    
    brain_dst = os.path.join(dst, os.path.basename(brain))
    if not os.path.exists(brain_dst): os.mkdir(brain_dst)
    
    for fl in fls_to_mv:
        shutil.move(fl, brain_dst)
