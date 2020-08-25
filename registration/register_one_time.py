#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:52:38 2019

@author: wanglab
"""

import os, sys
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.registration.register import elastix_command_line_call

# print(os.environ["SLURM_ARRAY_TASK_ID"])
# jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
jobid=0
imgs = ["/jukebox/LightSheetData/kocher-bee/volume_analysis/Grp16_2.575.tif",
        "/jukebox/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575.tif"]
mv = imgs[jobid]

fx = "/jukebox/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_rotate.tif"

out = mv[:-4]+"_elastix"
print(out)

mv = "/jukebox/LightSheetData/kocher-bee/volume_analysis/Grp16_2.575_elastix/result.0.tif"
out = "/jukebox/LightSheetData/kocher-bee/volume_analysis/Grp16_2.575_elastix_2"
if not os.path.exists(out): os.mkdir(out)

param_fld = "/jukebox/LightSheetData/kocher-bee/volume_analysis/parameter_files"
params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]

e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
