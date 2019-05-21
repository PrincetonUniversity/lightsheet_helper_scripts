#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 11:16:37 2019

@author: wanglab
"""

import tifffile

src = "/jukebox/wang/willmore/lightsheet/20190510_fiber_placement"

brains = ["/home/wanglab/mounts/wang/willmore/lightsheet/20190510_fiber_placement/m058_ventral_up/20190510_lw_m058_ventral_up_1d3x_488_555_017na_1hfds_z10um_100msec_resized_ch01_resampledforelastix.tif",
          "/home/wanglab/mounts/wang/willmore/lightsheet/20190510_fiber_placement/m058_dorsal_up/20190510_lw_m058_dorsal_up_1d3x_488_555_017na_1hfds_z10um_100msec_resized_ch01_resampledforelastix.tif",
          "/home/wanglab/mounts/wang/willmore/lightsheet/20190510_fiber_placement/m361_dorsal_up/20190510_lw_m361_dorsal_up_1d3x_488_555_017na_1hfds_z10um_100msec_resized_ch01_resampledforelastix.tif", 
          "/home/wanglab/mounts/wang/willmore/lightsheet/20190510_fiber_placement/m361_ventral_up/20190510_lw_m361_ventral_up_1d3x_488_555_017na_1hfds_z10um_100msec_resized_ch01_resampledforelastix.tif",
          "/home/wanglab/mounts/wang/willmore/lightsheet/20190510_fiber_placement/m363_dorsal_up/20190510_lw_m363_dorsal_up_1d3x_488_555_017na_1hfds_z10um_100msec_resized_ch01_resampledforelastix.tif",
          "/home/wanglab/mounts/wang/willmore/lightsheet/20190510_fiber_placement/m363_ventral_up/20190510_lw_m363_ventral_up_1d3x_488_555_017na_1hfds_z10um_100msec_resized_ch01_resampledforelastix.tif"]


for brain in brains:
    tifffile.imsave(brain, tifffile.imread(brain)[1:])
    