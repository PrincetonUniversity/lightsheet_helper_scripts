#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:52:38 2019

@author: wanglab
"""

import os, sys
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.registration.register import elastix_command_line_call

mv = "/jukebox/LightSheetData/brodyatlas/atlas/for_registration_to_lightsheet/WHS_SD_rat_T2star_v1.01_atlas.tif"

fx = "/jukebox/LightSheetData/rat-brodyprocessed/201910_tracing/z265/pbibawi_z265_ctbtracing_4x_647_017na_1hfds_z10um_150msec_20povlp_resized_ch00.tif"

out = "/jukebox/LightSheetData/rat-brodyprocessed/201910_tracing/z265/elastix"
if not os.path.exists(out): os.mkdir(out)

param_fld = "/jukebox/LightSheetData/brodyatlas/atlas/for_registration_to_lightsheet/rat_registration_parameter_folder"
params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]

e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
