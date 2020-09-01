#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 10:14:22 2020

@author: wanglab
"""

import os, sys
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.registration.register import elastix_command_line_call

if __name__ == "__main__":

    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    fxs = ["/jukebox/LightSheetData/kocher-bee/volume_analysis/Grp16_2.575.tif",
          "/jukebox/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575.tif"]
    fx = fxs[jobid]
    mv = "/jukebox/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575.tif"
    
    out = "/jukebox/LightSheetData/kocher-bee/volume_analysis/template_to_brain"
    out = os.path.join(out, os.path.basenmae(fx)[:-4]+"_elastix")
    print(out)
    
    if not os.path.exists(out): os.mkdir(out)
    
    param_fld = "/jukebox/LightSheetData/kocher-bee/volume_analysis/parameter_files"
    params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    
    e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)