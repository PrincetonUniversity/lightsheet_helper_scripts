#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:26:43 2019

@author: wanglab
"""

from tools.registration.register import elastix_command_line_call

fx = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"

mv = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"

out = "/home/wanglab/Desktop/aba_to_pma"


params = ["/jukebox/wang/pisano/Python/atlas/parameterfolder/Order1_Par0000affine.txt",
          "/jukebox/wang/pisano/Python/atlas/parameterfolder/Order2_Par0000bspline.txt"]

e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)