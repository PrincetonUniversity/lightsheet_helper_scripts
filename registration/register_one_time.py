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

fx = "/jukebox/wang/pisano/tracing_output/bl6_ts/20150804_tp_bl6_ts04/20150804_tp_bl6_ts04_488w_647_z3um_250msec_1hfds_resized_ch01_resampledforelastix.tif"
mv = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"

out = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125/20150804_tp_bl6_ts04/elastix_inverse_transform"
print(out)

if not os.path.exists(out): os.mkdir(out)

param_fld = "/jukebox/wang/zahra/python/BrainPipe/parameterfolder"
params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]

e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
