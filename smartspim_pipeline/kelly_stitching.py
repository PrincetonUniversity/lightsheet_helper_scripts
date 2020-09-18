#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:03:52 2020

@author: wanglab
"""

import os, subprocess

if __name__ == "__main__":
    
    #set paths
    src = "/jukebox/LightSheetData/wang-mouse/seagravesk"
    brains = [os.path.join(src,xx) for xx in os.listdir(src)]
    
    #paralellize across brains
    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    #voxel size
    xyz=(1.81,1.81,2)
    #select brain
    brain = brains[jobid]
    
    dst = os.path.join(brain, "Ex_785_Em_3", "stitched")
    if not os.path.exists(dst):
        print("\nStitching for brain %s" % os.path.basename(brain))
        #import
        print(subprocess.check_output(["terastitcher", "-1", "--volin=%s" % (os.path.join(brain, "Ex_785_Em_3")),
            "--ref1=x", "--ref2=y", "--ref3=z", "--vxl1=%0.2f" % xyz[0], 
            "--vxl2=%0.2f" % xyz[1], "--vxl3=%0.2f" % xyz[2], "--projout=xml_import"])) #command based on my previous run
        #compute
        print("\nDisplacement computation, this may take a while... \n\n")
        print(subprocess.check_output(["terastitcher", "--displcompute",
             "--projin=%s" % (os.path.join(brain, "Ex_785_Em_3", "xml_import.xml")),
             "--subvoldim=100", "--projout=xml_displcomp"]))
        #project, threshold, place tiles
        print("\nProjecting and thresholding... \n\n")
        print(subprocess.check_output(["terastitcher", "--displproj", 
             "--projin=%s" % (os.path.join(brain, "Ex_785_Em_3", "xml_displcomp.xml"))]))
        #thresholding
        print(subprocess.check_output(["terastitcher", "--displthres",
             "--projin=%s" % (os.path.join(brain, "Ex_785_Em_3", "xml_displproj.xml")),
             "--projout=%s" % (os.path.join(brain, "Ex_785_Em_3", "xml_displthres.xml")),
             "--threshold=0.7"]))
        #tiles
        print(subprocess.check_output(["terastitcher", "--placetiles",
             "--projin=%s" % (os.path.join(brain, "Ex_785_Em_3", "xml_displthres.xml")), 
             "--projout=%s" % (os.path.join(brain, "Ex_785_Em_3", "xml_placetiles.xml")),
             "--algorithm=MST"]))
        #merging
        os.mkdir(dst) #make stitched directory
        print("\nBlending and merging...\n\n")
        print(subprocess.check_output(["terastitcher", " --merge",
             "--projin=%s" % (os.path.join(brain, "Ex_785_Em_3", "xml_placetiles.xml")),
             "--volout=%s" % dst,
             "--imout_depth=16", "--resolutions=0"]))
    else:
        print("\nDoes not need stitching")

        
