#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 16:57:14 2021

@author: wanglab
"""

from scipy.io import loadmat
import tifffile as tif, numpy as np, os, matplotlib.pyplot as plt, sys
import shutil 

def transformix_command_line_call(src, dst, transformfile):
    '''Wrapper Function to call transformix using the commandline, this can be time consuming
    
    Inputs
    -------------------
    src = volume path for transformation
    dst = folder to save file
    transformfile = final transform file from elastix registration
    
    '''
    from subprocess import check_output
    print ('Running transformix, this can take some time....\n')
    #sp.call(['transformix', '-in', src, '-out', dst, '-tp', transformfile])
    call = 'transformix -in {} -out {} -tp {}'.format(src, dst, transformfile)
    print(check_output(call, shell=True))
    print('Past transformix command line Call')      
            
    return

def change_transform_parameter_initial_transform(fl, initialtrnsformpth):
    '''
    (InitialTransformParametersFileName "NoInitialTransform")
    initialtrnsformpth = 'NoInitialTransform' or 'pth/to/transform.0.txt'
    '''
    fl1 = fl[:-5]+'0000.txt'
    with open(fl, 'r') as f, open(fl1, 'w') as new:
            for line in f:
                new.write('(InitialTransformParametersFileName "{}")\n'.format(initialtrnsformpth)) if 'InitialTransformParametersFileName' in line else new.write(line)
    #overwrite original transform file
    shutil.move(fl1, fl)
    return

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

#path to mat file
mat = "/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_15_27_24_f37077_demonstrator_20171011/Ex_785_Em_3/corrected/sliding_diff_peak_find_99percentile_test20200806_all_coord_format2.mat"
pnts = loadmat(mat)["cell_centers_orig_coord"].astype(int)
#for resize dimensions
downsized = "/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_15_27_24_f37077_demonstrator_20171011/Ex_785_Em_3/reg_downsized_for_atlas.tif"
downsized = tif.imread(downsized) #sagittal
zd,yd,xd = downsized.shape #sagittal
#reorient pnts
pnts_sag = np.array([[xx[2],xx[1],xx[0]] for xx in pnts])
#get full size dims
stitched = "/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_15_27_24_f37077_demonstrator_20171011/Ex_785_Em_3/stitched/RES(7565x5726x2961)/098550/098550_104079"
y,z = tif.imread(os.path.join(stitched, os.listdir(stitched)[0])).shape #sagittal
x = len([xx for xx in os.listdir(stitched) if ".tif" in xx]) #sagittal
f = ((zd/z),(yd/y),(xd/x))
downsized_pnts_sag = np.array([[xx[0]*f[0],xx[1]*f[1],xx[2]*f[2]] for xx in pnts_sag]).astype(int)

#transform
#make into transformix-friendly text file
transformed_dst = "/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_15_27_24_f37077_demonstrator_20171011/Ex_785_Em_3/points"
if not os.path.exists(transformed_dst): os.mkdir(transformed_dst)
pretransform_text_file = create_text_file_for_elastix(downsized_pnts_sag, transformed_dst)
transformfiles = ["/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_15_27_24_f37077_demonstrator_20171011/elastix_inverse_transform/TransformParameters.0.txt",                  
                  "/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_15_27_24_f37077_demonstrator_20171011/elastix_inverse_transform/TransformParameters.1.txt",                  
                  "/jukebox/wang/seagravesk/lightsheet/201710_cfos_left_side_only_registration/f37077_demons/elastix_inverse_transform/cellch_f37077_demonstrator_20171011_790_015na_1hfsds_z5um_1000msec/f37077_demonstrator_20171011_488_015na_1hfsds_z5um_150msec_resized_ch00_resampledforelastix_atlas2reg/TransformParameters.0.txt",
                  "/jukebox/wang/seagravesk/lightsheet/201710_cfos_left_side_only_registration/f37077_demons/elastix_inverse_transform/cellch_f37077_demonstrator_20171011_790_015na_1hfsds_z5um_1000msec/f37077_demonstrator_20171011_488_015na_1hfsds_z5um_150msec_resized_ch00_resampledforelastix_atlas2reg/TransformParameters.1.txt"]
#copy over elastix files
transformfiles = modify_transform_files(transformfiles, transformed_dst) 
change_transform_parameter_initial_transform(transformfiles[0], 'NoInitialTransform')
#run transformix on points
points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
#convert registered points into structure counts
converted_points = unpack_pnts(points_file, transformed_dst)
#%%
atl_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
atl = tif.imread(atl_pth)
z,y,x = atl.shape
#check
if isinstance(converted_points, str):
    converted_points = np.load(converted_points)
arr=converted_points.astype(int)
cell=np.zeros((z,y,x)) #init cellmap
miss = 0
for pnt in arr:
    z,y,x=pnt
    try:
        cell[z,y,x] = 1
    except:
        miss+=1

plt.imshow(cell[300])