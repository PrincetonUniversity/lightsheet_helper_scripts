#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:52:38 2019

@author: wanglab
"""

import os
from tools.registration.register import elastix_command_line_call

mv = "/jukebox/LightSheetData/brodyatlas/atlas/for_registration_to_lightsheet/WHS_SD_rat_T2star_v1.01_atlas.tif"

fx = "/jukebox/LightSheetData/brodyatlas/atlas/2019_meta_atlas/median_image.tif"

out = "/jukebox/LightSheetData/brodyatlas/atlas/for_registration_to_lightsheet/waxholm_to_pra"
if not os.path.exists(out): os.mkdir(out)

param_fld = "/jukebox/LightSheetData/brodyatlas/atlas/for_registration_to_lightsheet/rat_registration_parameter_folder"
params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]

e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)