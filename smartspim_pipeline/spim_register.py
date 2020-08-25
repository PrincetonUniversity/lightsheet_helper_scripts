#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:01:17 2020

@author: wanglab
"""

import sys,os
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.registration.register import elastix_command_line_call

#takes 4 command line arguments max
print(sys.argv)
stepid = int(sys.argv[1])
src = str(sys.argv[2]) #folder to stitched images
reg = str(sys.argv[3]) #folder fo registration channel, e.g. Ex_488_Em_0
try:
    cell = str(sys.argv[4]) #folder for cell channel e.g. Ex_642_Em_2
except:
    cell = False

param_fld = "/jukebox/wang/zahra/python/BrainPipe/parameterfolder" #change if using rat
#registration vol to atlas
# cell = True #no separate cell/registration volume
# src = "/jukebox/LightSheetTransfer/tp/20200701_12_55_28_20170207_db_bl6_crii_rpv_01/"

if stepid == 0:
    reg = os.path.join(src, reg)
    mv = os.path.join(reg, "downsized_for_atlas.tif")
    print("\nPath to downsized vol for registration to atlas: %s" % mv)
    fx = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
    print("\nPath to atlas: %s" % fx)
    out = os.path.join(os.path.dirname(src), "elastix")
    if not os.path.exists(out): os.mkdir(out)
    
    params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    #run
    e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)

    if cell:
        #cell vol to registration vol
        print("\nCell channel specified: %s" % cell)
        mv = os.path.join(src, cell+"/downsized_for_atlas.tif")
        fx = os.path.join(src, reg+"/downsized_for_atlas.tif")
        
        out = os.path.join(src, "elastix/%s_to_%s" % (cell, reg))
        if not os.path.exists(out): os.mkdir(out)
        
        params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
        #run
        e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)

elif stepid == 1:
    #atlas to registration vol
    #inverse transform
    fx = os.path.join(src, reg+"/downsized_for_atlas.tif")
    mv = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
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
        mv = os.path.join(src, reg+"/downsized_for_atlas.tif")
        fx = os.path.join(src, cell+"/downsized_for_atlas.tif")
        
        out = os.path.join(src, "elastix_inverse_transform/%s_to_%s" % (reg, cell))
        if not os.path.exists(out): os.mkdir(out)
        
        params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
        #run
        e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
        
        