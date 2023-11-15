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
    stepid = int(sys.argv[1]) # 0 = normal registration, 1 = inverse registration
    src = str(sys.argv[2]) #folder to main image folder
    reg = str(sys.argv[3]) #folder fo registration channel, e.g. Ex_488_Em_0
    reg_vol = str(sys.argv[4])
    try:
        cell = str(sys.argv[5]) #folder for cell channel e.g. Ex_642_Em_2
        cell_vol = str(sys.argv[6])
    except:
        cell = False

    if cell == "NA":
        cell = False
    # try:
    #     species = str(sys.argv[5]) #species to know for registration parameters
    #     param_fld = "/jukebox/wang/ahoag/brainpipe/parameterfolder" #change if using rat
    # except:
    param_fld = "/jukebox/wang/ahoag/brainpipe/parameterfolder" #change if using rat
    

    atl_name = str(sys.argv[7])
    atl = ""
    if atl_name == "PMA":
        atl = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif" # PMA 20 micron isotropic
    elif atl_name == "Allen":
        atl = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif" # allen atlas 25 micron isotropic
    elif atl_name == "cb":
        #atl = "/jukebox/LightSheetTransfer/atlas/cb_sagittal_atlas_20um_iso.tif"
        atl = "/jukebox/LightSheetTransfer/pma_cropped_atlas.tif"
    elif atl_name == "PRA":
        atl = "/jukebox/brody/lightsheet/atlasdir/mPRA.tif"
        param_fld = "/jukebox/brody/lightsheet/elastix_params"
    elif atl_name == "cz":
        atl = "/jukebox/witten/Chris/data/clearmap2/utilities/allen-atlas-cz/average_template_25_sagittal_forDVscans_cz.tif"
    elif atl_name == "hem":
        atl = "/jukebox/LightSheetTransfer/hem_sagittal_atlas.tif"
    else:
        raise ValueError("Specified atlas does not exist")

    print(src)
    print(reg_vol)
    assert os.path.exists(param_fld)
    if stepid == 0:
        print("Doing normal registration")
        mv = os.path.join(reg, reg_vol)
        print(mv)
        assert os.path.exists(mv)
        print("\nPath to downsized vol for registration to atlas: %s" % mv)
        fx = atl
        assert os.path.exists(fx)
        print("\nPath to atlas: %s" % fx)
        out = os.path.join(src, "elastix_488_to_atl")
        if not os.path.exists(out): os.mkdir(out)
        params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
        #run
        e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)

        if cell:
            #cell vol to registration vol
            print("\nCell channel specified: %s" % cell)
            mv = os.path.join(cell, cell_vol)
            fx = os.path.join(reg, reg_vol)          

            out = os.path.join(src, "elastix/cell_to_reg")
            if not os.path.exists(out): os.mkdir(out)
            
            params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
            #run
            e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)

    elif stepid == 1:
        print("Doing inverse registration")
        #atlas to registration vol
        #inverse transform
        fx = os.path.join(reg, reg_vol)
        mv = atl
        assert os.path.exists(fx)
        assert os.path.exists(mv)
        print("\nPath to downsized vol for inverse registration to atlas: %s" % fx)
        print("\nPath to atlas: %s" % mv)
        out = os.path.join(src, "elastix_inverse_transform")
        if not os.path.exists(out): os.mkdir(out)
        
        params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
        #run
        e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
        
        #registration vol to cell vol
        #inverse transform
        if cell:
            print("\nCell channel specified: %s" % cell)
            mv = os.path.join(reg, reg_vol)
            fx = os.path.join(cell, cell_vol)
            
            assert os.path.exists(fx)
            assert os.path.exists(mv)
            
            out = os.path.join(src, "elastix_inverse_transform", "reg_to_cell")
            if not os.path.exists(out): os.mkdir(out)
            
            params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
            #run
            e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
        
        
