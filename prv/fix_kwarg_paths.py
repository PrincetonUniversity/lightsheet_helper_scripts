#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 17:18:30 2019

@author: wanglab
"""

import os
from tools.utils.io import load_kwargs, save_kwargs

fld = "/jukebox/wang/pisano/tracing_output/retro_4x/"

brains = ["20180313_jg_bl6f_prv_23",
         "20180322_jg_bl6f_prv_28",
         "20180313_jg_bl6f_prv_20",
         "20180323_jg_bl6f_prv_31",
         "20180322_jg_bl6f_prv_27",
         "20180313_jg_bl6f_prv_21",
         "20180326_jg_bl6f_prv_34",
         "20180326_jg_bl6f_prv_36",
         "20180323_jg_bl6f_prv_30",
         "20180326_jg_bl6f_prv_32",
         "20180322_jg_bl6f_prv_29",
         "20180326_jg_bl6f_prv_33",
         "20180326_jg_bl6f_prv_37",
         "20180322_jg_bl6f_prv_26",
         "20180313_jg_bl6f_prv_25",
         "20180326_jg_bl6f_prv_35",
         "20180313_jg_bl6f_prv_24"]

for brain in brains:
    print(brain)
    src = os.path.join(fld, brain)
    kwargs = load_kwargs(src, update_dict = True, systemdirectory = "/jukebox/")
    save_kwargs(**kwargs)

    