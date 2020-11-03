#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 10:14:22 2020

@author: wanglab
"""

import os, sys
from subprocess import check_output
import subprocess as sp


def elastix_command_line_call(fx, mv, out, parameters, fx_mask=False):
    """Wrapper Function to call elastix using the commandline, this can be time consuming
    
    Inputs
    -------------------
    fx = fixed path (usually Atlas for "normal" noninverse transforms)
    mv = moving path (usually volume to register for "normal" noninverse transforms)
    out = folder to save file
    parameters = list of paths to parameter files IN ORDER THEY SHOULD BE APPLIED
    fx_mask= (optional) mask path if desired
    
    Outputs
    --------------
    ElastixResultFile = ".tif" or ".mhd" result file
    TransformParameterFile = file storing transform parameters
    
    """
    e_params=["elastix", "-f", fx, "-m", mv, "-out", out]
    if fx_mask: e_params=["elastix", "-f", fx, "-m", mv, "-fMask", fx_mask, "-out", out]
    
    ###adding elastix parameter files to command line call
    for x in range(len(parameters)):
        e_params.append("-p")
        e_params.append(parameters[x])
    print("Elastix Command:\n{}\n...".format(e_params))    
    
    #set paths
    TransformParameterFile = os.path.join(out, "TransformParameters.{}.txt".format((len(parameters)-1)))
    ElastixResultFile = os.path.join(out, "result.{}.tif".format((len(parameters)-1)))
    
    #run elastix: 
    try:                
        print("Running Elastix, this can take some time....\n")
        sp.call(e_params)#sp_call(e_params)#
        print(out,"Past Elastix Commandline Call")
    except RuntimeError as e:
        print(out,"\n***RUNTIME ERROR***: {} Elastix has failed. Most likely the two images are too dissimiliar.\n".format(e.message))
        pass      
    if os.path.exists(ElastixResultFile) == True:    
        print(out,"Elastix Registration Successfully Completed\n")
    #check to see if it was MHD instead
    elif os.path.exists(os.path.join(out, "result.{}.mhd".format((len(parameters)-1)))) == True:    
        ElastixResultFile = os.path.join(out, "result.{}.mhd".format((len(parameters)-1)))
        print(out,"Elastix Registration Successfully Completed\n")
    else:
        print(out, "\n***ERROR***Cannot find elastix result file\n: {}".format(ElastixResultFile))
        return
        
    return ElastixResultFile, TransformParameterFile

def transformix_command_line_call(src, dst, transformfile):
    """Wrapper Function to call transformix using the commandline, this can be time consuming
    
    Inputs
    -------------------
    src = volume path for transformation
    dst = folder to save file
    transformfile = final transform file from elastix registration
    
    """
    
    print ("Running transformix, this can take some time....\n")
    #sp.call(["transformix", "-in", src, "-out", dst, "-tp", transformfile])
    call = "transformix -in {} -out {} -tp {}".format(src, dst, transformfile)
    print(check_output(call, shell=True))
    print("Past transformix command line Call")      
            
    return

if __name__ == "__main__":
    #run
    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    #set paths
    src = "/jukebox/LightSheetData/kocher-bee/volume_analysis/volumes_downsized_to_template"
    dst = "/jukebox/LightSheetData/kocher-bee/volume_analysis/"
    #brains
    brs = ["D07ret_2.575.tif", "C40iso_2.575.tif", "D18grp_2.575.tif",
       "C09ret_2.575.tif", "A01iso_2.575step.tif", "D47ret_2.575.tif",
       "C18ret_2.575.tif", "C30iso_2.575.tif", "B04ret_2.575.tif",
       "C03grp_2.575step.tif", "D05iso_2.575.tif", "D16grp_2.575.tif",
       "B02ret_2.575.tif", "D14grp_2.575.tif", "C28iso_2.576.tif",
       "D21grp_2.575.d2.tif", "D08ret_2.575.tif", "C07iso_2.575.tif",
       "C19ret_2.575.tif", "C37grp_2.575.tif", "C04grp_2.575.tif",
       "D42ret_2.575.tif", "C16ret_2.575.tif", "B01ret_2.575step.tif",
       "D12grp_2.575.tif", "B03ret_2.575.tif", "C08iso_2.575.tif",
       "B07iso_2.575.tif", "C11ret_2.575.tif", "C25grp_2.575.tif",
       "C13ret_2.575.tif", "D01grp_2.575.tif", "D24iso_2.575.tif",
       "D04grp_2.575.tif", "C05grp_2.575.tif", "D41ret_2.575.tif",
       "C12ret_2.575.tif", "B08iso_2.575.tif", "C15ret_2.575.tif",
       "D17grp_2.575.tif", "D40ret_2.575.tif", "C14ret_2.575.tif",
       "D27grp_2.575.tif", "B12iso_2.575.tif", "C33grp_2.575.tif",
       "B11iso_2.575.tif", "C29iso_2.575.tif", "B05ret_2.575.tif"]
    #'fixed' imagess
    fxs = [os.path.join(src,xx) for xx in brs]
    #array job
    fx = fxs[jobid]
    ##template brain
    mv = os.path.join(dst,"template/Bombus45_2.575umstep_rotate_croppedZ.tif")
    #output dir
    out = os.path.join(dst,"template_to_brain")
    out = os.path.join(out, os.path.basename(fx)[:-4]+"_elastix")
    print(out)
    #make sure directory tree exist
    if not os.path.exists(os.path.dirname(out)): os.mkdir(os.path.dirname(out))
    if not os.path.exists(out): os.mkdir(out)
    #set parameter files
    param_fld = os.path.join(dst,"parameter_files")
    params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
    #run
    e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)

    #transform atlas
    transform = transformfiles
    #change transform file output datatype
    with open(transform, "r") as file:
        filedata = file.read()
        # Replace the target string
        #make sure outtput is float32
        filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
        #for atlas transforms
        filedata = filedata.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)')
        # Write the file out again
        with open(transform, "w") as file:
          file.write(filedata)
    #transform annotation file to experimental brain      
    ann = "/jukebox/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment_croppedZ.tif"
    dst = out
    #run
    transformix_command_line_call(ann, dst, transform)