#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 13:00:36 2020

@author: wanglab
"""


import os, shutil

src = "/jukebox/wang/mkislin/lightsheet_brains/201903_cntnap2_tsc1_ai148/todo"
dst = src
#format of raw data filename you want to convert this to
temp = "10-50-20_UltraII_raw_RawDataStack[00 x 00]_C00_xyz-Table %s_UltraII Filter000%s.ome.tif"
brains = ["201707_mk61_488_647_014na_1hfsds_z10um_300msec_ch00",
            "201707_mk61_488_647_014na_1hfsds_z10um_300msec_ch01",
            "201707_mk62_488_647_014na_1hfsds_z10um_300msec_ch00",
            "201707_mk62_488_647_014na_1hfsds_z10um_300msec_ch01",
            "201707_mk63_488_647_014na_1hfsds_z10um_300msec_ch00",
            "201707_mk63_488_647_014na_1hfsds_z10um_300msec_ch01"]

for brain in brains:
    print("\n"+brain+"\n")
    #set brain dst
    brain_dst = os.path.join(dst, brain)
    #join source paths
    brain = os.path.join(src, brain)
    fls = [os.path.join(brain, xx) for xx in os.listdir(brain)]; fls.sort()
    for jjj, fl in enumerate(fls):
        zpln = fl[-9:-4]
        renm = temp % (zpln,0) #rename to fit raw data
        shutil.copy(fl, os.path.join(brain_dst, renm))
        if jjj%50==0: print("\n**********plane %s copied!*********\n" % jjj)
