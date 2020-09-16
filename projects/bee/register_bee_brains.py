#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 10:14:22 2020

@author: wanglab
"""

import os, sys
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.registration.register import elastix_command_line_call, transformix_command_line_call

if __name__ == "__main__":

    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    src = "/jukebox/LightSheetData/kocher-bee/volume_analysis"
    fxs = [os.path.join(src,"Grp16_2.575.tif"),
          os.path.join(src,"IsoYellow_2.575.tif"),
          os.path.join(src,"retc17_2.575umstep.tif")]
    fx = fxs[jobid]
    mv = os.path.join(src,"template/Bombus45_2.575umstep_rotate_croppedZ.tif")
    
    out = os.path.join(src,"template_to_brain")
    out = os.path.join(out, os.path.basename(fx)[:-4]+"_elastix")
    print(out)
    
    if not os.path.exists(os.path.dirname(out)): os.mkdir(os.path.dirname(out))
    if not os.path.exists(out): os.mkdir(out)
    
    param_fld = os.path.join(src,"parameter_files")
    params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)

    #transform atlas
    transform = transformfiles[-1]
    with open(transform, "r") as file:
        filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
        filedata = filedata.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)')
        # Write the file out again
        with open(transform, "w") as file:
          file.write(filedata)
          
    ann = "/jukebox/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment_croppedZ.tif"
    dst = out
    transformix_command_line_call(ann, dst, transform)


    #test
    # fx = os.path.join(src,"IsoYellow_2.575.tif")
    # mv = os.path.join(src,"Grp16_2.575.tif")
    # out = "/home/wanglab/Desktop/beetest"
    # if not os.path.exists(out): os.mkdir(out)
    # param_fld = os.path.join(src,"parameter_files")
    # params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    # e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)