#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 16:57:14 2021

@author: wanglab
"""

from scipy.io import loadmat
import tifffile as tif, numpy as np, os, matplotlib.pyplot as plt, sys, pandas as pd
import shutil, itertools 
from collections import Counter

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

def fast_scandir(dirname):
    """ gets all folders recursively """
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders

def transformed_pnts_to_allen_helper_func(arr, ann, order = "XYZ"):
    """Function to transform given array of indices and return the atlas pixel ID from the annotation file
    
    Input
    --------------
    numpy array of Nx3 dimensions corresponding to ***XYZ*** coordinates
    ann = numpy array of annotation file
    order = "XYZ" or "ZYX" representing the dimension order of arr"s input
    
    Returns
    -------------
    Pixel value at those indices of the annotation file
    """        
    ########procecss
    pnt_lst=[]; badpntlst = []
    for i in range(len(arr)):
        try:        
            pnt=[int(x) for x in arr[i]]
            if order == "XYZ": pnt_lst.append(ann[pnt[2], pnt[1], pnt[0]]) ###find pixel id; arr=XYZ; ann=ZYX
            elif order == "ZYX": pnt_lst.append(ann[pnt[0], pnt[1], pnt[2]]) ###find pixel id; arr=ZYX; ann=ZYX
        except IndexError:
            badpntlst.append([pnt[2], pnt[1], pnt[0]]) #ZYX
            pass ######THIS NEEDS TO BE CHECKED BUT I BELIEVE INDEXES WILL BE OUT OF 
    sys.stdout.write("\n*************{} points do not map to atlas*********\n".format(len(badpntlst))); sys.stdout.flush()
    return pnt_lst

def count_structure_lister(allen_id_table, *args):
    """Function that generates a pd table of structures where contour detection has been observed
    Inputs:
        allen_id_table=annotation file as np.array
        *args=list of allen ID pixel values ZYX
    """
    #make dictionary of pixel id:#num of the id
    cnt = Counter()
    for i in args:
        cnt[i]+=1
    #generate df + empty column
    if type(allen_id_table) == str: allen_id_table = pd.read_excel(allen_id_table, index_col=None) #df=allen_id_table.assign(count= [0]*len(allen_id_table)) #add count columns to df
    df=allen_id_table
    df["cell_count"]=0
    #populate cell count in dataframe
    for pix_id, count in cnt.items():
        df.loc[df.id==pix_id, "cell_count"]=count
    return df

if __name__ == "__main__":
    
    #pipeline to loop thru all brains with detected cells
    src = "/jukebox/LightSheetData/wang-mouse/seagravesk"
    #list of animals
    # animals = ["20200916_19_25_35_f37080_mouse1_20171015"]
    animals = ["20200810_13_10_58_f37080_mouse2_20171015_slow_focus", "20200901_14_20_11_f37073_mouse1_20171010", "20200901_15_27_24_f37077_demonstrator_20171011",
                "20200901_16_34_36_m37112_observer_20171010", "20200901_17_43_19_f37106_mouse2_20171011", "20200902_13_04_58_m37083_demonstrator_20171018",
                "20200902_14_13_54_f37078_observer_20171014","20200902_16_50_10_f37073_mouse2_20171010", "20200902_17_53_21_m37111_demonstrator_20171012",
                "20200903_11_51_58_m37111_observer_20171012", "20200903_12_51_14_m37110_demonstrator_20171016", "20200903_14_01_08_m37072_demonstrator_20171008",
                "20200909_13_16_53_m37109_mouse2_20171018", "20200909_14_54_45_f37104_observer_20171016", "20200909_15_59_15_f37105_observer_20171012",
                "20200909_17_01_50_m37112_demonstrator_20171010","20200911_14_01_24_m37079_mouse2_20171014","20200911_14_57_36_m37081_observer_20171014",
                "20200911_16_45_49_m37081_demonstrator_20171014","20200911_17_46_51_f37105_demonstrator_20171012","20200916_14_37_19_f37107_demonstrator_20171007",
                "20200916_16_42_48_m37109_mouse1_20171018","20200916_19_25_35_f37080_mouse1_20171015",
                "20200917_12_27_13_m37072_observer_20171008","20200917_13_24_45_m37071_demonstrator_20171006","20200917_14_22_17_f37104_demonstrator_20171016",
                "20200921_12_14_34_m37110_observer_20171016","20200921_13_13_19_f37070_demonstrator_20171007","20200921_14_31_11_m37071_observer_20171006",
                "20200921_15_37_17_f37107_observer_20171007","20200921_16_45_26_m37113_mouse1_20171007","20200923_10_56_13_f37070_observer_20171007",
                "20200923_12_05_45_f37077_observer_20171011","20200923_13_21_40_m37083_observer_20171018","20200923_14_37_00_m37079_mouse1_20171014",
                "20200924_13_33_09_m37113_mouse2_20171007","20200924_14_33_12_f37106_mouse1_20171011","20201102_16_29_12_m37110_demonstrator_20171016"]
    missing = []
    for animal in animals:
        try:
            name = animal[18:]
            print("\n\n **********{}********* \n\n".format(name))
            #path to mat file
            try: #for differences in folder structure
                mat = os.path.join(src, animal, "Ex_785_Em_3/corrected/sliding_diff_peak_find_99percentile_test20200806_all_coord_format2.mat")
                pnts = loadmat(mat)["cell_centers_orig_coord"].astype(int)
            except:
                correctedpth = os.path.join(src, animal, "Ex_785_Em_3/corrected")
                matpth = fast_scandir(correctedpth)[-1]
                mat = os.path.join(matpth, "sliding_diff_peak_find_99percentile_test20200806_all_coord_format2.mat")
                if not os.path.exists(mat):
                    #move 1 folder up, for brain with extra folder
                    try:
                        mat = os.path.join(fast_scandir(correctedpth)[-2], "sliding_diff_peak_find_99percentile_test20200806_all_coord_format2.mat")
                    except:
                        mat = os.path.join(fast_scandir(correctedpth)[-3], "sliding_diff_peak_find_99percentile_test20200806_all_coord_format2.mat")
                pnts = loadmat(mat)["cell_centers_orig_coord"].astype(int)
            #for resize dimensions
            downsized = os.path.join(src, animal, "Ex_785_Em_3/reg_downsized_for_atlas.tif")
            downsized = tif.imread(downsized) #sagittal
            zd,yd,xd = downsized.shape #sagittal
            #reorient pnts
            pnts_sag = np.array([[xx[2],xx[1],xx[0]] for xx in pnts])
            #get full size dims
            stitched = os.path.join(src, animal, "Ex_785_Em_3/stitched")
            stitched = fast_scandir(stitched)[-1] #gets to the end of directory tree
            y,z = tif.imread(os.path.join(stitched, os.listdir(stitched)[0])).shape #sagittal
            x = len([xx for xx in os.listdir(stitched) if ".tif" in xx]) #sagittal
            f = ((zd/z),(yd/y),(xd/x))
            downsized_pnts_sag = np.array([[xx[0]*f[0],xx[1]*f[1],xx[2]*f[2]] for xx in pnts_sag]).astype(int)
            
            #transform
            #make into transformix-friendly text file
            wanglabserver = "/jukebox/wang/seagravesk/lightsheet/smartspim_points_202107"
            transformed_dst = os.path.join(wanglabserver, "{}_points".format(name))
            if not os.path.exists(transformed_dst): os.mkdir(transformed_dst)
            pretransform_text_file = create_text_file_for_elastix(downsized_pnts_sag, transformed_dst)
            #get old atlas transform fils
            oldr = "/jukebox/wang/seagravesk/lightsheet/201710_cfos_left_side_only_registration"
            oldinvertransform = os.path.join(oldr, "{}/elastix_inverse_transform".format(name[:13]))
            cellchtransform = [os.path.join(oldinvertransform,xx) for xx in os.listdir(oldinvertransform) if "cellch" in xx][0]
            atlas2reg = [os.path.join(cellchtransform, xx) for xx in os.listdir(cellchtransform) if "atlas2reg" in xx][0]
            trlst = [os.path.join(atlas2reg,xx) for xx in os.listdir(atlas2reg) if "TransformParameters" in xx]; trlst.sort()
            transformfiles = [os.path.join(src, animal, "elastix_inverse_transform/TransformParameters.0.txt"),                  
                              os.path.join(src, animal, "elastix_inverse_transform/TransformParameters.1.txt"),
                              trlst[0],trlst[1]]
            #copy over elastix files
            transformfiles = modify_transform_files(transformfiles, transformed_dst) 
            change_transform_parameter_initial_transform(transformfiles[0], 'NoInitialTransform')
            #run transformix on points
            points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
            #convert registered points into structure counts
            converted_points = unpack_pnts(points_file, transformed_dst)
            #align to annotation
            ann_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif"
            id_table = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table.xlsx"
            point_lst = transformed_pnts_to_allen_helper_func(np.load(converted_points), 
                        tif.imread(ann_pth), order = "ZYX")
            #zmd added 20190312 - these should be in order of points inputted from raw space
            np.save(os.path.join(transformed_dst, "annotation_pixel_value_coordinates.npy"), point_lst)            
            df = count_structure_lister(id_table, *point_lst).fillna(0)
            df = df.drop(columns = ["Unnamed: 0"])
            df.to_csv(os.path.join(transformed_dst, os.path.basename(id_table).replace(".xlsx", "")+"_with_anatomical_assignment_of_cell_counts.csv"),
                      index = None)
            with open(os.path.join(transformed_dst, "info_from_zahras_code.txt"), "w") as fl:
                fl.write("FILES AND PATHS USED TO GENERATE ANALYSIS:\n")
                fl.write("\n")
                fl.write("path_to_csv_file: {}\n".format(os.path.join(transformed_dst, os.path.basename(id_table).replace(".xlsx", "")+"_with_anatomical_assignment_of_cell_counts.csv")))
                fl.write("ann_pth: {}\n".format(ann_pth))
                fl.write("id_table: {}\n".format(id_table))
                fl.write("dst: {}\n".format(transformed_dst))
                fl.write("cell_centers_path: {}\n".format(mat))
                fl.write("\n")
                fl.write("\n")
                fl.write("Date code was run: {}\n".format(datetime.datetime.today().strftime("%Y-%m-%d")))
                fl.close()
        except:
            missing.append(animal)
# #%%
# #check
# atl_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
# atl = tif.imread(atl_pth)
# z,y,x = atl.shape
# #check
# if isinstance(converted_points, str):
#     converted_points = np.load(converted_points)
# arr=converted_points.astype(int)
# cell=np.zeros((z,y,x)) #init cellmap
# miss = 0
# for pnt in arr:
#     z,y,x=pnt
#     try:
#         cell[z,y,x] = 1
#     except:
#         miss+=1

# plt.imshow(cell[300])