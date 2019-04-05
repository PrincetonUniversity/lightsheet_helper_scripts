#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 11:29:49 2019

@author: wanglab
"""

from skimage.external import tifffile as tif
import numpy as np
import pandas as pd
import multiprocessing as mp
from collections import Counter

def transformed_pnts_to_allen_helper_func(arr, ann, order = 'XYZ'):
    '''Function to transform given array of indices and return the atlas pixel ID from the annotation file
    
    Input
    --------------
    numpy array of Nx3 dimensions corresponding to ***XYZ*** coordinates
    ann = numpy array of annotation file
    order = 'XYZ' or 'ZYX' representing the dimension order of arr's input
    
    Returns
    -------------
    Pixel value at those indices of the annotation file, maintains order if NO BAD MAPPING
    '''        
    ########procecss
    pnt_lst = []
    arr = list(arr)
    for i in range(len(arr)):
        pnt = [x for x in arr[i]]
        if order == 'XYZ': pnt_lst.append(ann[pnt[2], pnt[1], pnt[0]]) ###find pixel id; arr=XYZ; ann=ZYX
        elif order == 'ZYX': pnt_lst.append(ann[pnt[0], pnt[1], pnt[2]]) ###find pixel id; arr=ZYX; ann=ZYX
    return pnt_lst
        
def count_structure_lister(allen_id_table, args):
    '''Function that generates a pd table of structures 
    Inputs:
        allen_id_table = annotation file
        args = dictionary of of allen ID pixel values ZYX
    '''
    #set dataframe
    #generate df + empty column
    if type(allen_id_table) == str: allen_id_table = pd.read_excel(allen_id_table)
    df = allen_id_table
    df['voxels_in_structure'] = 0
        
    #make dictionary of pixel id:#num of the id
    for arg in args:
        arg = args[arg]
        cnt = Counter()
        for i in arg:
            cnt[i] += 1
        
        #populate cell count in dataframe
        for pix_id, count in cnt.items():
            df.loc[df.id == pix_id, 'voxels_in_structure'] = count

    return df

#%%

if __name__ == "__main__":

    ann = '/home/wanglab/mounts/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif'
    id_table = '/home/wanglab/mounts/LightSheetTransfer/atlas/allen_atlas/allen_id_table.xlsx'
    
    ann = tif.imread(ann)
    
    structs = np.unique(ann)
    
    areas = {}
    
    #heavy
    for struct in structs:
        areas[struct] = np.where(ann == struct)
        
    ##test
    #nz = areas[1085.0]    
    #pos = transformed_pnts_to_allen_helper_func(np.asarray(zip(*[nz[0], nz[1], nz[2]])), ann, "ZYX")    
    #
    #tdf = count_structure_lister(id_table, *pos)
    
    #heavy
    pos = {}
    for area in areas:
        nz = areas[area]
        pos[area] = transformed_pnts_to_allen_helper_func(zip(nz[0], nz[1], nz[2]), ann, "ZYX")
        
    
    df = count_structure_lister(id_table, pos)
    
    df.to_csv('/home/wanglab/Desktop/allen_id_table_w_voxel_counts.csv')
    df = df.drop(columns = ["cell_count", "Unnamed: 0"])
    df.to_excel('/home/wanglab/Desktop/allen_id_table_w_voxel_counts.xlsx')