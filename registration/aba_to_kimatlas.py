#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:26:43 2019

@author: wanglab
"""

from tools.registration.register import elastix_command_line_call
import os

#aba --> kimatlas a.k.a Paxinos/Allen hybrid
#in terms of the transform
mv = "/jukebox/LightSheetTransfer/atlas/kim_atlas/KimRef_tissue_volume_4brainpipe.tif"

fx = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"

out = "/jukebox/wang/ahoag/aba_to_kimatlas"
if not os.path.exists(out): os.mkdir(out)

params = ["/jukebox/wang/pisano/Python/atlas/parameterfolder/Order1_Par0000affine.txt",
          "/jukebox/wang/pisano/Python/atlas/parameterfolder/Order2_Par0000bspline.txt"]

e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)