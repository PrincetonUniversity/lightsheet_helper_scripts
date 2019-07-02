#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:47:05 2019

@author: wanglab
"""

import os, re

pth = "/home/wanglab/mounts/LightSheetData/mallarino/ricardo/20190702/raw_data/190702_20190313_edu_171_20190610_mallarino_1d3x_488_008na_z25um_1hfds_100ms_40povlp_sagittal_11-03-29"

imgs = [os.path.join(pth, xx) for xx in os.listdir(pth) if "UltraII_raw_RawDataStack" in xx]; imgs.sort()

regexpression = r'(.*)(?P<y>\d{2})(.*)(?P<x>\d{2})(.*C+)(?P<ls>[0-9]{1,2})(.*Z+)(?P<z>[0-9]{1,4})(.*r)(?P<ch>[0-9]{1,4} )(.ome.tif)'

reg = re.compile(regexpression)
matches = list(map(reg.match, imgs)) #matches[0].groups()

##find index of z,y,x,ch in a match str
z_indx = matches[0].span('z')
try:
    y_indx = matches[0].span('y')
    x_indx = matches[0].span('x')
    tiling = True
except IndexError:
    y_indx = 1
    x_indx = 1
    tiling = False

chs=[]; [chs.append(matches[i].group('ch')[:]) for i in range(len(matches)) if matches[i].group('ch')[:] not in chs]
zplns=[]; [zplns.append(matches[i].group('z')) for i in range(len(matches)) if matches[i].group('z') not in zplns]; zplns.sort()
