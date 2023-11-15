#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:01:17 2020

@author: wanglab
"""

import sys,os
sys.path.append("/jukebox/wang/sanjeev/BrainPipe")
from tools.registration.register import elastix_command_line_call

if __name__ == '__main__':
    #takes 6 command line arguments max
    print(sys.argv)
    src = str(sys.argv[1]) #folder to main image folder
    dst = str(sys.argv[2]) #folder fo registration channel, e.g. Ex_488_Em_0
    # try:
    #     species = str(sys.argv[5]) #species to know for registration parameters
    #     param_fld = "/jukebox/wang/ahoag/brainpipe/parameterfolder" #change if using rat
    # except:
    param_fld = "/jukebox/wang/ahoag/brainpipe/parameterfolder" #change if using rat
    # atl = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif" # PMA 20 micron isotropic
    atl = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif" # allen atlas 25 micron isotropic
    print(src)
    print(dst)
    assert os.path.exists(param_fld)

    print("Doing normal registration")
    mv = os.path.join(src, "downsized_for_atlas.tif")
    print(mv)
    assert os.path.exists(mv)
    print("\nPath to downsized vol for registration to atlas: %s" % mv)
    fx = os.path.join(dst, "downsized_for_atlas.tif")
    assert os.path.exists(fx)
    print("\nPath to atlas: %s" % fx)
    out = os.path.join(os.path.dirname(src), "elastix")
    if not os.path.exists(out): os.mkdir(out)
    params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    #run
    e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)

        
        
