#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 18:33:02 2019

@author: wanglab
"""

import os, numpy as np, time, sys
from collections import Counter
from skimage.external import tifffile
import matplotlib.pyplot as plt
import matplotlib.colors
from io import BytesIO
import zipfile
import pandas as pd

def transformed_pnts_to_allen_helper_func(arr, ann, order = "XYZ"):
    """Function to transform given array of indices and return the atlas pixel ID from the annotation file
    
    Input
    --------------
    numpy array of Nx3 dimensions corresponding to ***XYZ*** coordinates
    ann = numpy array of annotation file
    order = "XYZ" or "ZYX" representing the dimension order of arr"s input
    
    Returns
    -------------
    Pixel value at those indices of the annotation file, maintains order if NO BAD MAPPING
    """        
    ########procecss
    pnt_lst = []; badpntlst = []
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
    if type(allen_id_table) == str: allen_id_table = pd.read_excel(allen_id_table) #df=allen_id_table.assign(count= [0]*len(allen_id_table)) #add count columns to df
    df=allen_id_table
    df["cell_count"]=0
    
    #populate cell count in dataframe
    for pix_id, count in cnt.iteritems():
        df.loc[df.id==pix_id, "cell_count"]=count

    return df

def read_roi(fileobj):
    """
    points = read_roi(fileobj)

    Read ImageJ"s ROI format. Points are returned in a nx2 array. Each row
    is in [row, column] -- that is, (y,x) -- order.
    """
    # This is based on:
    # http://rsbweb.nih.gov/ij/developer/source/ij/io/RoiDecoder.java.html
    # http://rsbweb.nih.gov/ij/developer/source/ij/io/RoiEncoder.java.html

    SPLINE_FIT = 1
    DOUBLE_HEADED = 2
    OUTLINE = 4
    OVERLAY_LABELS = 8
    OVERLAY_NAMES = 16
    OVERLAY_BACKGROUNDS = 32
    OVERLAY_BOLD = 64
    SUB_PIXEL_RESOLUTION = 128
    DRAW_OFFSET = 256

    class RoiType:
        POLYGON = 0
        RECT = 1
        OVAL = 2
        LINE = 3
        FREELINE = 4
        POLYLINE = 5
        NOROI = 6
        FREEHAND = 7
        TRACED = 8
        ANGLE = 9
        POINT = 10

    def get8():
        s = fileobj.read(1)
        if not s:
            raise IOError("readroi: Unexpected EOF")
        return ord(s)

    def get16():
        b0 = get8()
        b1 = get8()
        return (b0 << 8) | b1

    def get32():
        s0 = get16()
        s1 = get16()
        return (s0 << 16) | s1

    def getfloat():
        v = np.int32(get32())
        return v.view(np.float32)

    magic = fileobj.read(4)
    if magic != b"Iout":
        raise ValueError("Magic number not found")
    version = get16()

    # It seems that the roi type field occupies 2 Bytes, but only one is used
    roi_type = get8()
    # Discard second Byte:
    get8()

    if roi_type not in [RoiType.FREEHAND, RoiType.POLYGON, RoiType.RECT, RoiType.POINT]:
        raise NotImplementedError("roireader: ROI type %s not supported" % roi_type)

    top = get16()
    left = get16()
    bottom = get16()
    right = get16()
    n_coordinates = get16()
    x1 = getfloat()
    y1 = getfloat()
    x2 = getfloat()
    y2 = getfloat()
    stroke_width = get16()
    shape_roi_size = get32()
    stroke_color = get32()
    fill_color = get32()
    subtype = get16()
    if subtype != 0:
        raise NotImplementedError("roireader: ROI subtype %s not supported (!= 0)" % subtype)
    options = get16()
    arrow_style = get8()
    arrow_head_size = get8()
    rect_arc_size = get16()
    position = get32()
    header2offset = get32()

    if roi_type == RoiType.RECT:
        if options & SUB_PIXEL_RESOLUTION:
            return np.array(
                [[y1, x1], [y1, x1+x2], [y1+y2, x1+x2], [y1+y2, x1]],
                dtype=np.float32)
        else:
            return np.array(
                [[top, left], [top, right], [bottom, right], [bottom, left]],
                dtype=np.int16)

    if options & SUB_PIXEL_RESOLUTION:
        getc = getfloat
        points = np.empty((n_coordinates, 2), dtype=np.float32)
        fileobj.seek(4*n_coordinates, 1)
    else:
        getc = get16
        points = np.empty((n_coordinates, 2), dtype=np.int16)

    points[:, 1] = [getc() for i in range(n_coordinates)]
    points[:, 0] = [getc() for i in range(n_coordinates)]

    if options & SUB_PIXEL_RESOLUTION == 0:
        points[:, 1] += left
        points[:, 0] += top

    return points

def read_roi_zip(fname, include_roi_name=False, verbose=True):
    """Wrapper for reading zip files generated from ImageJ (FIJI)
    
    include_roi_name (optional) 
        if true: returns list of (roi_name, contour)
        roi_name=z,y,x
        useful for extracting z (NOTE: ImageJ has one-based numerics vs Python w zero-based numerics)
    """        
    try:
        if not include_roi_name:
            with zipfile.ZipFile(fname) as zf:
                return [read_roi(zf.open(n)) for n in zf.namelist()]
                                                    
        if include_roi_name:
            with zipfile.ZipFile(fname) as zf:
                return [(n, read_roi(zf.open(n))) for n in zf.namelist()]
    
    #hack to try and keep 
    except ValueError:
        import sys
        sys.stdout.write("possibly error with the ROI file but probably ok")

#%%                                     
if __name__ == "__main__":
    
    fld = "/jukebox/LightSheetData/witten-mouse/201903_Alex_NeuroPixels_CMDiI/20190314_ibl_cm_dii_01"
    #load in ROIS - clicked in horizontal volume
    roi_pth = "[PATH TO ROI ZIP FILE]"
    
    start = time.time()
      
    brain = os.path.basename(fld)
    
    #check if path exists obv
    if os.path.exists(roi_pth):
            
        zyx_rois = np.asarray([[int(yy) for yy in xx.replace(".roi", "").split("-")] for xx in read_roi_zip(roi_pth, include_roi_name=True)])
            
        #make merged image
        zyx = np.asarray([str((int(xx[0]), int(xx[1]), int(xx[2]))) for xx in zyx_rois])
        
        #make destination path
        dst = os.path.join(fld, "ROI_merged_to_atlas")
        if not os.path.exists(dst): os.mkdir(dst)
        
        #export coordinates
        if os.path.exists(os.path.join(dst, "{}_allen_coordinates.txt".format(brain))): os.remove(os.path.join(dst, "{}_Allen_coordinates.txt".format(brain)))
        with open(os.path.join(dst, "{}_allen_coordinates.txt".format(brain)), "a") as txt:
            txt.write("Allen Atlas CCF coordinates (in zyx):\n%s" % zyx_rois)
        
        #atlas (horizontal)
        atl = tifffile.imread("/jukebox/LightSheetData/witten-mouse/atlas/average_template_25_sagittal_forDVscans.tif")
        atl = np.transpose(atl, [2, 1, 0])
        atl_cnn = np.zeros_like(atl)
        
        for i in range(zyx_rois.shape[0]):
            atl_cnn[zyx_rois[i][0], zyx_rois[i][1], zyx_rois[i][2]] = 1
     
        merged = np.stack([atl, atl_cnn, np.zeros_like(atl)], -1)
        tifffile.imsave(os.path.join(dst, "{}_ROI_merged_to_Allen_horizontal.tif".format(brain)), merged)
        
        coronal = np.transpose(merged, [1, 0, 2, 3]) #make coronal sections
        
        #save out coronal sections - based on the fact that you click 5 points in each site in HORIZONTAL sectionss
        tifffile.imsave(os.path.join(dst, "{}_ROI_merged_to_Allen_coronal.tif".format(brain)), coronal)
        
        #doing a max projection, in case you just want to look at that
        maxip = np.max(coronal, 0)
        
        alpha = 0.6 #determines transparency, don"t need to alter
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #color to visualise
        plt.imshow(maxip[...,0], "gist_yarg")
        plt.imshow(maxip[...,1], cmap, alpha = 0.6)
        plt.axis("off")
        plt.savefig(os.path.join(dst, "{}_ROI_merged_to_Allen_coronal_maxip.pdf".format(brain)), dpi = 300)
        
        print("\n\ntook {} seconds to make merged maps for {}\n".format(time.time()-start, brain))
        
        #make allen structure LUT
        #go from horiztonal to sag
        zyx_rois = np.asarray([[xx[2], xx[1], xx[0]] for xx in zyx_rois])
        
        #convert to structure
        annotation_file = "/jukebox/LightSheetData/witten-mouse/atlas/annotation_template_25_sagittal_forDVscans.tif"
        ann = tifffile.imread(annotation_file)
        points = transformed_pnts_to_allen_helper_func(list(zyx_rois), ann, order = "ZYX")    
        
        #make dataframe
        lut_path = "/jukebox/LightSheetData/witten-mouse/atlas/allen_id_table_w_voxel_counts.xlsx"
        df = count_structure_lister(lut_path, *points)
        df.to_excel(os.path.join(dst, "{}_Allen_structures.xlsx".format(brain)))