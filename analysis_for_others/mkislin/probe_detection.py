#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:27:56 2019

@author: wanglab
"""

import tifffile

impth = "/home/wanglab/mounts/wang/oostland/lightsheet/m26/elastix/marlies_m26_1d3x_488_555_008na_1hfds_z5um_150msec_resized_ch01/result.tif"

im = tifffile.imread(impth)