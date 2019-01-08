#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 11:16:03 2019

@author: wanglab
"""
import os, sys
sys.path.append('/jukebox/wang/zahra/clearmap_cluster_copy/')
from ClearMap.Analysis.Label import countPointsInRegions

import numpy as np, pandas as pd, xlrd

def make_table_of_transformed_cells(src):
    ''' 
    incorporates new atlas and look up table for anatomical analysis of cfos data done using:
        https://github.com/PrincetonUniversity/clearmap_cluster
    NOTE: this assumes the clearmap transform has run correctly, and uses transformed cells
    '''
    
    #first, check if cells where transformed properly
    #TRANSFORMED cell dataframe
    points = np.load(os.path.join(src, "clearmap_cluster_output/cells_transformed_to_Atlas.npy"))
    raw = np.load(os.path.join(src, "clearmap_cluster_output/cells.npy"))
    
    print(points.shape, raw.shape)
    
    if points.shape == raw.shape:
        intensities = np.load(os.path.join(src, "clearmap_cluster_output/intensities.npy"))
        
        #LUT
        ann = '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif'
        ann_lut = '/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx'
        
        #open LUT excel sheet
        wb = xlrd.open_workbook(ann_lut)
        lut = wb.sheet_by_index(0)
        
        #Table generation:
        ##################
        #With integrated weigths from the intensity file (here raw intensity):
        ids, counts = countPointsInRegions(points, labeledImage = ann, intensities = intensities, intensityRow = 0);
        #mapping
        #NOTE: REMOVES 'basic cell groups and region since LUT doesn't have a value for that. FIX!??!
        id2name = dict((row[3].value, row[1]) for row in (lut.row(r) for r in xrange(lut.nrows)))
       
        table = {}
        table["id"] = ids[1:]
        table["counts"] = counts[1:].astype('int64')
        
        table["name"] = [id2name[i_d].value for i_d in ids[1:]]
        
        pd.DataFrame.from_dict(table, orient = "columns").to_csv(os.path.join(src, 'Annotated_counts_intensities.csv'))
            
        #Without weigths (pure cell number):
        #mapping
        #NOTE: REMOVES 'basic cell groups and region since LUT doesn't have a value for that. FIX!??!
        id2name = dict((row[3].value, row[1]) for row in (lut.row(r) for r in xrange(lut.nrows)))
       
        table = {}
        table["id"] = ids[1:]
        table["counts"] = counts[1:].astype('int64')
        
        table["name"] = [id2name[i_d].value for i_d in ids[1:]]
        
        pd.DataFrame.from_dict(table, orient = "columns").to_csv(os.path.join(src, 'Annotated_counts.csv'))
            
        print ('Analysis Completed') 
    else:
        print ('Cells not transformed properly, check transformix (aka step 6) and try again')
    
#%%
        
if __name__ == '__main__':
    #goal is to transform cooridnates, voxelize based on number of cells and overlay with reigstered cell signal channel...
    #inputs
#    src = '/jukebox/wang/Jess/lightsheet_output/201810_cfos/processed/dadult_pc_crus1_3'
    pth = '/home/wanglab/mounts/wang/Jess/lightsheet_output/201810_cfos/processed'
    for src in os.listdir(pth):
        make_table_of_transformed_cells(os.path.join(pth, src))
