#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:52:38 2019

@author: wanglab
"""

import os
from tools.registration.register import elastix_command_line_call

mv = "/home/wanglab/mounts/wang/pisano/tracing_output/antero/20170308_tp_bl6F_cri_2x_03/registration_to_aba/20170308_tp_bl6F_cri_2x_03_488_555_647_0005na_1hfsds_z3um_225msec_resized_ch01_resampledforelastix.tif"

fx = "/home/wanglab/mounts/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"

out = "/home/wanglab/mounts/wang/zahra/h129_contra_vs_ipsi/rtn/20170308_tp_bl6f_cri_2x_03_reg"
if not os.path.exists(out): os.mkdir(out)


params = ["/home/wanglab/mounts/wang/zahra/python/lightsheet_py3/parameterfolder/Order1_Par0000affine.txt",
          "/home/wanglab/mounts/wang/zahra/python/lightsheet_py3/parameterfolder/Order2_Par0000bspline.txt"]

e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)