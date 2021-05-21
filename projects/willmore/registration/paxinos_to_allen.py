#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:01:17 2020

@author: wanglab
"""

import sys,os
sys.path.append("/jukebox/wang/ahoag/brainpipe")
from tools.registration.register import elastix_command_line_call

if __name__ == '__main__':
    #takes 6 command line arguments max
    param_fld = "/jukebox/wang/ahoag/brainpipe/parameterfolder" #change if using rat
    # atl = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif" # PMA 20 micron isotropic
    paxinos_atl = "/jukebox/LightSheetTransfer/atlas/kim_atlas/KimRef_tissue_volume_4brainpipe.tif" # Paxinos atlas 10 micron x 100 micron x 10 microns
    allen_atl = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif" # allen atlas 25 micron isotropic

    assert os.path.exists(param_fld)
    assert os.path.exists(paxinos_atl)
    assert os.path.exists(allen_atl)

    # Lindsay wants to transform paxinos points to allen space so 
    # we want the inverse transformation with allen as the moving image
    # and paxinos as the fixed image 
    fx = paxinos_atl
    mv = allen_atl
    print("Doing inverse registration")
    #atlas to registration vol
    #inverse transform

    print("\nPath to Fixed image: %s" % fx)
    print("\nPath to Moving image: %s" % mv)
    out = os.path.join("/jukebox/LightSheetTransfer/atlas/kim_atlas",
        "elastix_inverse_transform_allenatlas")
    if not os.path.exists(out): os.mkdir(out)
    
    params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    #run
    e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
