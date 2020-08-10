#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 15:39:28 2020

@author: wanglab
"""


import os, sys, numpy as np, shutil
from scipy.io import loadmat, savemat
import pickle

def transform_points(src, dst, transformfiles, resample_points=False):
    """
    
    Inputs
    ---------
    src = numpy file consiting of nx3 (ZYX points)
    dst = folder location to write points
    transformfiles = 
        list of all elastix transform files used, and in order of the original transform****
    resample_points = [original_dims, resample_dims] if there was resampling done, use this here
    param_dictionary_for_reorientation = param_dictionary for lightsheet package to use for reorientation
    """
    #load
    cells = np.load(src) if src[-3:] == "npy" else np.array([loadmat(src)["all_z_coord_together"][0],
                                                                  loadmat(src)["all_y_coord_together"][0],
                                                                  loadmat(src)["all_x_coord_together"][0]]).T
    
    #reorient - test to make sure this works
    cells[:,[0,2]] = cells[:,[2,0]] #horizontal to saggittal 
            
    #optionally resample points
    if resample_points:
        original_dims, resample_dims = resample_points
        cells = points_resample(cells, original_dims, resample_dims)
    
    #generate text file
    pretransform_text_file = create_text_file_for_elastix(cells, dst)
    
    #copy over elastix files
    transformfiles = modify_transform_files(transformfiles, dst)
   
    #run transformix on points
    points_file = point_transformix(pretransform_text_file, transformfiles[-1], dst)
    
    #convert registered points into structure counts
    unpack_pnts(points_file, dst)   
    
    return
    
    
def create_text_file_for_elastix(src, dst):
    """
    
    Inputs
    ---------
    src = numpy file consiting of nx3 (ZYX points)
    dst = folder location to write points
    """
    
    print("This function assumes ZYX centers...")
    
    #setup folder
    if not os.path.exists(dst): os.mkdir(dst)
                                 
    #create txt file, with elastix header, then populate points
    pth=os.path.join(dst, "zyx_points_pretransform.txt")
    
    #load
    if type(src) == np.ndarray:
        arr = src
    else:
        arr = np.load(src) if src[-3:] == "npy" else loadmat(src)["cell_centers_orig_coord"]
    
    #convert
    stringtowrite = "\n".join(["\n".join(["{} {} {}".format(i[2], i[1], i[0])]) for i in arr]) ####this step converts from zyx to xyz*****
    
    #write file
    sys.stdout.write("writing centers to transfomix input points text file..."); sys.stdout.flush()
    with open(pth, "w+") as fl:
        fl.write("index\n{}\n".format(len(arr)))    
        fl.write(stringtowrite)
        fl.close()
    sys.stdout.write("...done writing centers\n"); sys.stdout.flush()
        
    return pth

def modify_transform_files(transformfiles, dst):
    """Function to copy over transform files, modify paths in case registration was done on the cluster, and tether them together
    
        Inputs
    ---------
    transformfiles = 
        list of all elastix transform files used, and in order of the original transform****
    
    """
    
    #new
    ntransformfiles = [os.path.join(dst, "order{}_{}".format(i,os.path.basename(xx))) for i,xx in enumerate(transformfiles)]
    
    #copy files over
    [shutil.copy(xx, ntransformfiles[i]) for i,xx in enumerate(transformfiles)]
    
    #modify each with the path
    for i,pth in enumerate(ntransformfiles):
        
        #skip first
        if i!=0:
            
            #read
            with open(pth, "r") as fl:
                lines = fl.readlines()
                fl.close()
            
            #copy
            nlines = lines
            
            #iterate
            for ii, line in enumerate(lines):
                if "(InitialTransformParametersFileName" in line:
                    nlines[ii] = "(InitialTransformParametersFileName {})\n".format(ntransformfiles[i-1])
            
            #write
            with open(pth, "w") as fl:
                for nline in lines:
                    fl.write(str(nline))
                fl.close()
        
    return ntransformfiles
    
    

def point_transformix(pretransform_text_file, transformfile, dst):
    """apply elastix transform to points
    
    
    Inputs
    -------------
    pretransform_text_file = list of points that already have resizing transform
    transformfile = elastix transform file
    dst = folder
    
    Returns
    ---------------
    trnsfrm_out_file = pth to file containing post transformix points
    
    """
    sys.stdout.write("\n***********Starting Transformix***********")
    from subprocess import check_output
    #set paths    
    trnsfrm_out_file = os.path.join(dst, "outputpoints.txt")
    
    #run transformix point transform
    call = "transformix -def {} -out {} -tp {}".format(pretransform_text_file, dst, transformfile)
    print(check_output(call, shell=True))
    sys.stdout.write("\n   Transformix File Generated: {}".format(trnsfrm_out_file)); sys.stdout.flush()
    return trnsfrm_out_file


def unpack_pnts(points_file, dst):
    """
    function to take elastix point transform file and return anatomical locations of those points
    
    Here elastix uses the xyz convention rather than the zyx numpy convention
    
    Inputs
    -----------
    points_file = post_transformed file, XYZ
    
    Returns
    -----------
    dst_fl = path to numpy array, ZYX
    
    """   

    #####inputs 
    assert type(points_file)==str
    point_or_index = 'OutputPoint'
    
    #get points
    with open(points_file, "r") as f:                
        lines=f.readlines()
        f.close()

    #####populate post-transformed array of contour centers
    sys.stdout.write("\n\n{} points detected\n\n".format(len(lines)))
    arr=np.empty((len(lines), 3))    
    for i in range(len(lines)):        
        arr[i,...]=lines[i].split()[lines[i].split().index(point_or_index)+3:lines[i].split().index(point_or_index)+6] #x,y,z
            
    #optional save out of points
    dst_fl = os.path.join(dst, "posttransformed_zyx_voxels.npy")
    np.save(dst_fl, np.asarray([(z,y,x) for x,y,z in arr]))
    
    #check to see if any points where found
    print("output array shape {}".format(arr.shape))
        
    return dst_fl

def points_resample(src, original_dims, resample_dims, verbose=False):
    """Function to adjust points given resizing by generating a transform matrix
    
    ***Assumes ZYX and that any orientation changes have already been done.***
    
    src: numpy array or list of np arrays of dims nx3
    original_dims (tuple)
    resample_dims (tuple)
    """
    src = np.asarray(src)
    assert src.shape[-1] == 3, "src must be a nx3 array"
    
    #initialize
    d1,d2=src.shape
    nx4centers=np.ones((d1,d2+1))
    nx4centers[:,:-1]=src
    
    #acount for resampling by creating transformmatrix
    zr, yr, xr = resample_dims
    z, y, x = original_dims
    
    #apply scale diff
    trnsfrmmatrix=np.identity(4)*(zr/float(z), yr/float(y), xr/float(x), 1)
    if verbose: sys.stdout.write("trnsfrmmatrix:\n{}\n".format(trnsfrmmatrix))
    
    #nx4 * 4x4 to give transform
    trnsfmdpnts=nx4centers.dot(trnsfrmmatrix) ##z,y,x
    if verbose: sys.stdout.write("first three transformed pnts:\n{}\n".format(trnsfmdpnts[0:3]))

    return trnsfmdpnts


def load_dictionary(pth):
    """simple function to load dictionary given a pth
    """
    kwargs = {};
    with open(pth, "rb") as pckl:
        kwargs.update(pickle.load(pckl))
        pckl.close()

    return kwargs

#%%
if __name__ == "__main__":
    
    
    ###NOTE CHECK TO ENSURE ACCOUNTING FOR INPUT RESAMPLING, and ORIENTATION CHANGE*****
    #inputs
    #numpy file consiting of nx3 (ZYX points) or if .mat file structure where zyx is called "cell_centers_orig_coord"
    src = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_785_Em_3/stitched/RES(7604x5720x3553)/102090/102090_120640/sliding_diff_peak_find_95percentile_test20200806_all_coord.mat"
    dst = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_785_Em_3/cells_95percentile_test20200806" # folder location to write points
    if not os.path.exists(dst): os.mkdir(dst)
    
    #when marking centers in the  "raw" full sized cfos channel. This will transform those centers into "atlas" space (in this case the moving image)
    #list of all elastix transform files used, and in order of the original transform****
    transformfiles = ["/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/elastix_inverse_transform/TransformParameters.0.txt", # this is auto = fixed image; atlas = moving image
                      "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/elastix_inverse_transform/TransformParameters.1.txt", # this is auto = fixed image; atlas = moving image
                      "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/elastix_inverse_transform/reg_to_sig785/TransformParameters.0.txt",  # this is cfos = fixed image; auto = moving image
                      "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/elastix_inverse_transform/reg_to_sig785/TransformParameters.1.txt"]  # this is cfos = fixed image; auto = moving image

    #optional resampling between fullsized and input to elastix
    original_dims = (5773,7574,3535) #sagittal 
    resample_dims = (638,739,448) #sagittal 
    resample_points = [original_dims, resample_dims]
    
    #apply
    transform_points(src, dst, transformfiles, resample_points)

#%%
    #check if mapping is correct, don't need to run if you are confident the script works
    import matplotlib.pyplot as plt, tifffile as tif
    cellspth = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_785_Em_3/cells_95percentile_test20200806/posttransformed_zyx_voxels.npy"
    cells = np.load(cellspth).astype(int)
    atlpth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
    # atlpth = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_785_Em_3/downsized_for_atlas.tif"
    atl = tif.imread(atlpth)
    #map all cells to atl
    cellvol = np.zeros_like(atl)
    # cellrs = cells.T[:3].T.astype(int)
    for i,cell in enumerate(cells):
        #only map if non negative coordinates
        if sum(cell>0)==3:
            if i%100000==0: print("%s cells mapped" % i)
            try:
                cellvol[cell[0],cell[1],cell[2]] = 255
            except:
                print("Cell coordinate out of bounds")
            
    plt.imshow(np.max(cellvol[250:300],axis=0))

#%%

    # Code for converting .npy to .mat
    f = "/some/path/.npy" # full path to the .npy that you want converted
    new_path = "some_other_path.mat" # full path to where the .mat file should be saved and what it should be named
    data = np.load(f) # Load the .npy matrix into the workspace
    save_dict = dict(x=data) # convert it to a dictionary
    savemat(new_path, save_dict) # save the dictionary as a .mat file
    