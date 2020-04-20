#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:35:14 2020

@author: wanglab
"""

import os, shutil

src = "/jukebox/LightSheetData/rat-brody/processed/201910_tracing/"
dst = "/jukebox/LightSheetTransfer/brody"
#format of raw data filename you want to convert this to
temp = "10-50-16_UltraII_raw_RawDataStack[00 x 00]_C00_xyz-Table %s_UltraII Filter000%s.ome.tif"
brains = ["z265", "z266", "z267", "z268", "z269"]

for brain in brains:
    print("\n"+brain+"\n")
    #set brain dst
    brain_dst = os.path.join(dst, brain)
    #join source paths
    brain = os.path.join(src, brain)
    if not os.path.exists(brain_dst): os.mkdir(brain_dst)
    #find fullsizedatafld
    fszfld = os.path.join(brain, "full_sizedatafld")
    #get all the channels stitched
    chs = [xx for xx in os.listdir(fszfld) if "msec" in xx]; chs.sort()
    for iii, ch  in enumerate(chs): #iii here becomes channel number
        print("\n**********channel %s*********\n" % iii)
        chfld = os.path.join(fszfld, ch)
        fls = [os.path.join(chfld, xx) for xx in os.listdir(chfld)]; fls.sort()
        for jjj, fl in enumerate(fls):
            zpln = fl[-9:-4]
            renm = temp % (zpln,iii) #rename to fit raw data
            shutil.copy(fl, os.path.join(brain_dst, renm))
            if jjj%10==0: print("\n**********plane %s copied!*********\n" % jjj)
