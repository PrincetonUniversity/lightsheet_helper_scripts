#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:19:28 2019

@author: wanglab
"""

import os, numpy as np, sys
from concurrent.futures import ProcessPoolExecutor
from skimage.external import tifffile
sys.path.append("/jukebox/wang/ahoag/brainpipe")
from tools.utils.io import makedir
from tools.registration.register import change_interpolation_order, transformix_command_line_call
from tools.registration.transform_list_of_points import modify_transform_files
from scipy.ndimage.interpolation import zoom

def process_slice(z):
    """ A function that is run in parallel below. 
    This function zooms out a single z plane to raw resolution
    and saves the zoomed out z plane. This function can read the 
    global variables from the if __name__ section so those 
    do not need to be passed to the function.
    """
    zplane = bigdvann[z]
    savename = os.path.join(save_dst, "annotation_Z%04d.tif" % z)
    if os.path.exists(savename):
        print("Already z plane # {}".format(z))
        return
    arr = zoom(zplane, (ap0/float(ap1), ml0/float(ml1)), order=0)
    tifffile.imsave(savename, arr, compress = 6)
    print("\nmade z plane # {}".format(z))
    return "success"

if __name__ == '__main__':
    # Read in the path to the brain whose raw data you want the atlas aligned to
    step = sys.argv[1]
    raw_dir = sys.argv[2]
    elastix_atlas_to_auto_dir = sys.argv[3]
    elastix_auto_to_cell_dir = sys.argv[4]
    output_dir = sys.argv[5]

    print(f"raw_dir is: {raw_dir}")
    print(f"raw_dir is: {raw_dir}")
    assert os.path.exists(raw_dir)
    assert len(os.listdir(raw_dir)) > 0
    #set which annotation atlas you want to align to raw space. e.g. Allen or PMA
    ann = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels.tif" # Princeton Mouse Annotation Atlas, 16 bit

    # Set up output directories    
    makedir(output_dir) # Makes it if it doesn't exist already

    aldst = os.path.join(output_dir, "transformed_annotations")
    makedir(aldst)

    save_dst = os.path.join(aldst, "single_tifs")
    makedir(save_dst)
    print(f"Output directory structure all set")
    sys.stdout.flush()

    if step == 'step1': # copy files and run transformix
        # Paths to the transform parameters from running elastix
        # Paths must be to the INVERSE transform files, i.e. the opposite direction to what you would logically think
        # reg -> cell, transform 0
        a2r0 = os.path.join(elastix_auto_to_cell_dir,"TransformParameters.0.txt")
        # reg -> cell, transform 1
        a2r1 = os.path.join(elastix_auto_to_cell_dir,"TransformParameters.1.txt")
        # atlas -> reg, transform 0
        r2s0 = os.path.join(elastix_atlas_to_auto_dir,"TransformParameters.0.txt")
        # atlas -> reg, transform 1
        r2s1 = os.path.join(elastix_atlas_to_auto_dir,"TransformParameters.1.txt")

        ## transformix
        # first copy over files to new location
        transformfiles = modify_transform_files(transformfiles=[a2r0, a2r1, r2s0, r2s1], dst = aldst)
        # change order of interpolation in new copied files
        [change_interpolation_order(xx,0) for xx in transformfiles]

        #changes the parameter in the transform files that outputs 16bit images instead
        for fl in transformfiles:# Read in the file
            with open(fl, "r") as file:
                filedata = file.read()
            # Replace the target string
            filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
            # Write the file out again
            with open(fl, "w") as file:
              file.write(filedata)
        print("Copied over transform files and modified them")

        # run transformix - makes the transformed atlas as result.tif in aldst
        print("running transformix")
        transformed_ann = os.path.join(aldst, "result.tif")
        if not os.path.exists(transformed_ann):
            transformix_command_line_call(ann, aldst, transformfiles[-1])
            print("ran transformix successfully!")
        else:
            print("transformix was already run before")
        sys.stdout.flush()

    if step == 'step2':
        #now zoom out - this is heavy!
        print("Now zooming out in z first. This can take a while")
        transformed_ann = os.path.join(aldst, "result.tif")
        tann = tifffile.imread(transformed_ann)
        pl0 = tifffile.imread(os.path.join(raw_dir, os.listdir(raw_dir)[0])) # Z=0 plane
        dv0,ap0,ml0 = len(os.listdir(raw_dir)), pl0.shape[0], pl0.shape[1] # dimensions of original raw data, z, y, z
        ml1,ap1,dv1 = tann.shape # dimensions of downsized raw data
        #scale in dv only first and rotate to hor orientation
        # bigdvann = np.flip(np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2),0) #ALSO REVERSE IN Z FOR VENTRAL TO DORSAL IMAGING
        bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)
        #remove np.flip if imaging_request_1g D --> V
        print("zoomed out in DV first")
        sys.stdout.flush()

        #now rotate and scale each in ap and ml
        to_process = list(range(0,len(bigdvann)))
        print(to_process)
        print("Rotating and scalling each zplane in AP in ML")
        with ProcessPoolExecutor(max_workers=8) as executor:
            for job in executor.map(process_slice,to_process):
                try:
                    print(job)
                except Exception as exc:
                    print(f'generated an exception: {exc}')
            

