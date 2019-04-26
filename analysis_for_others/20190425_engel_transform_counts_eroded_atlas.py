#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:57:37 2019

@author: wanglab
"""

import os
import numpy as np, pandas as pd, tifffile 

def labelPoints(points, labeledImage, level = None, collapse = None):
   """ borrowed/modified/cleaned up from Clearmap """
   
   x = points[:,2]; #ASSUMES TRANSFORMED CELLS ARE Z,Y,X??!?!??!
   y = points[:,1];
   z = points[:,0]; 
    
   nPoint = x.size;    
    
   pointLabels = np.zeros(nPoint, 'int32');
    
   labelImage = tifffile.imread(labeledImage)
   dsize = labelImage.shape;
    
   for i in range(nPoint):
       
       if x[i] >= 0 and x[i] < dsize[0] and y[i] >= 0 and y[i] < dsize[1] and z[i] >= 0 and z[i] < dsize[2]:
            pointLabels[i] = labelImage[int(x[i]), int(y[i]), int(z[i])];
            
   return pointLabels
            
def countPointsInRegions(points, labeledImage, intensities = None, intensityRow = 0, level= None, sort = True, 
                         returnIds = True, returnCounts = False, collapse = None):    
    """ borrowed/modified/cleaned up from Clearmap """
    
    pointLabels = labelPoints(points, labeledImage, level = level, collapse = collapse); 
    
    if intensities is None:
        ll, cc = np.unique(pointLabels, return_counts = True);
        cci = None;
    else:
        if intensities.ndim > 1:
            intensities = intensities[:,intensityRow];
   
        ll, ii, cc = np.unique(pointLabels, return_counts = True, return_inverse = True);
        cci = np.zeros(ll.shape);
        for i in range(ii.shape[0]):
             cci[ii[i]] += intensities[i]
    
    #cc = numpy.vstack((ll,cc)).T;
    if sort:
        ii = np.argsort(ll);
        cc = cc[ii];
        ll = ll[ii];
        if not cci is None:
            cci = cci[ii];

    if returnIds:
        if cci is None:
            return ll, cc
        else:
            if returnCounts:
                return ll, cc, cci;
            else:
                return ll, cci
    else:
        if cci is None:
            return cc;
        else:
            if returnCounts:
                return cc, cci;
            else:
                return cci;
            
def make_table_of_transformed_cells(src):
    """ 
    incorporates new atlas and look up table for anatomical analysis of cfos data done using:
        https://github.com/PrincetonUniversity/clearmap_cluster
    NOTE: this assumes the clearmap transform has run correctly, and uses transformed cells
    """
    
    print(src)
    #first, check if cells where transformed properly
    #TRANSFORMED cell dataframe
    try:
        points = np.load(os.path.join(src, "clearmap_cluster_output/cells_transformed_to_Atlas.npy"))
        raw = np.load(os.path.join(src, "clearmap_cluster_output/cells.npy"))
        
        print(points.shape, raw.shape)
        
        if points.shape == raw.shape:
            intensities = np.load(os.path.join(src, "clearmap_cluster_output/intensities.npy"))
            
            #LUT
            ann = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_75um_erosion.tif"
            ann_lut = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/ARA2_annotation_info_w_voxel_counts.csv"
            
            #open LUT csv
            lut = pd.read_csv(ann_lut)
            
            #Table generation:
            ##################
            #With integrated weigths from the intensity file (here raw intensity):
            ids, intensities = countPointsInRegions(points, labeledImage = ann, intensities = intensities, intensityRow = 0)
            #keep them together to modify together later
            ids_intensities = list(zip(ids.astype("float64"), intensities))
            
            #mapping
            id2name = {}; id2acronym = {}; id2parentacr = {}; id2voxcount = {}
            #NOTE: REMOVES "basic cell groups and region since LUT doesn"t have a value for that. FIX!??!
            for iid in lut.id.values:
                id2name[lut.loc[(lut.id == iid), "id"].values[0]] = lut.loc[(lut.id == iid), "name"].values[0]
                id2parentacr[lut.loc[(lut.id == iid), "id"].values[0]] = lut.loc[(lut.id == iid), "parent_acronym"].values[0]
                id2acronym[lut.loc[(lut.id == iid), "id"].values[0]] = lut.loc[(lut.id == iid), "acronym"].values[0]
                id2voxcount[lut.loc[(lut.id == iid), "id"].values[0]] = lut.loc[(lut.id == iid), "voxels_in_structure"].values[0]
            
            lut_ids = lut.id.values[1:]
           
            table = {}
            table["id"] = [i_d[0] for i_d in ids_intensities if i_d[0] in lut_ids]
            table["intensity"] = [i_d[1] for i_d in ids_intensities if i_d[0] in lut_ids]
            
            #dropping structures that do not map to LUT
            table["name"] = [id2name[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            table["acronym"] = [id2acronym[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            table["parent_acronym"] = [id2parentacr[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            table["voxels_in_structure"] = [id2voxcount[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            
            pd.DataFrame.from_dict(table, orient = "columns").to_csv(os.path.join(src, "Annotated_counts_intensities_75um_erosion.csv"), index = None)
                
            #Without weigths (pure cell number):
            ids, counts = countPointsInRegions(points, labeledImage = ann, intensities = None)
            ids_counts = list(zip(ids.astype("float64"), counts))
            
            table = {}
            table["id"] = [i_d[0] for i_d in ids_counts if i_d[0] in lut_ids]
            table["counts"] = [i_d[1] for i_d in ids_counts if i_d[0] in lut_ids]
            
            #dropping structures that do not map to LUT
            table["name"] = [id2name[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            table["acronym"] = [id2acronym[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            table["parent_acronym"] = [id2parentacr[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            table["voxels_in_structure"] = [id2voxcount[i_d] for i_d in ids[1:] if i_d in id2name.keys()]
            
            pd.DataFrame.from_dict(table, orient = "columns").to_csv(os.path.join(src, "Annotated_counts_75um_erosion.csv"), index = None)
                
            print ("\n Analysis Completed\n") 
        else:
            print ("\n Cells not transformed properly, check transformix (aka step 6) and try again\n")
    except:
        print("\n Path for transformed cells doesn't exist\n")
    
#%%
        
if __name__ == "__main__":
    #goal is to transform cooridnates, voxelize based on number of cells and overlay with reigstered cell signal channel...
    #inputs
#    src = 
#    make_table_of_transformed_cells(src)

    pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/processed"
    
    for src in os.listdir(pth):
        make_table_of_transformed_cells(os.path.join(pth, src))