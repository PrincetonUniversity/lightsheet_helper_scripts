#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:49:11 2019

@author: wanglab
"""
from tools.utils.io import load_kwargs
import tifffile
import matplotlib.pyplot as plt
pth = "/home/wanglab/mounts/wang/pisano/tracing_output/aav/20171130_pcp2_dcn_02"

kwargs = load_kwargs(pth)

cellvol = [xx for xx in kwargs["volumes"] if xx.ch_type == "injch" and xx.channel == "00"][0]

regpth = cellvol.ch_to_reg_to_atlas

img = tifffile.imread(regpth)

plt.imshow(img[300]*10)
