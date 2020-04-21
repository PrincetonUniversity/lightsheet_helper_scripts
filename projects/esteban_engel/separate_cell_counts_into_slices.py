#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 17:17:29 2019

@author: wanglab
"""

import os
import numpy as np, pandas as pd, tifffile, json, seaborn as sns

    # Now write the function to get all progeny of an input nodename
def get_progeny(dic,parent_structure,progeny_list):
    """ 
    ---PURPOSE---
    Get a list of all progeny of a structure name.
    This is a recursive function which is why progeny_list is an
    argument and is not returned.
    ---INPUT---
    dic                  A dictionary representing the JSON file 
                         which contains the ontology of interest
    parent_structure     The structure
    progeny_list         The list to which this function will 
                         append the progeny structures. 
    """
    if 'msg' in list(dic.keys()): dic = dic['msg'][0]
    
    name = dic.get('name')
    children = dic.get('children')
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get('name')
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get('name')
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 


def generate_cell_density_per_slice_df(src, dst, ann, structures, ids_all_structures, zr_per_struct):
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
        
        #NOTE THAT POINTS HERE ARE SAVED AS XYZ - making them sagittal
        points = np.asarray([(int(xx[2]), int(xx[1]), int(xx[0])) for xx in points])
        raw = np.load(os.path.join(src, "clearmap_cluster_output/cells.npy"))        
        print(points.shape, raw.shape)
        
        if points.shape == raw.shape:
            
            #so if we consider our sub volumes to comprimise 20 planes each...we need to find the number of the voxels that
            #structure occupies in those 20 planes and divide the total cell count for that structure in those 20 planes
            #to get the density value
            zr_per_slice = [np.arange(zr[0], zr[1]+19, 20).astype(int) for zr in zr_per_struct] #get the ranges for the slices per structure
            
            counts_per_structure = []
            for strt, zr in enumerate(zr_per_slice): #per structure
                counts_per_slice = []
                for i,slc in enumerate(zr): #for each slice in the structure
                    if i < len(zr)-1: #so that the 2nd to last number in the array is considered as the beginning of the last slice
                        counts = 0 #init counter
                        counts = [pnt for pnt in points if pnt[0] >= slc and
                          pnt[0] < zr[i+1] and pnt[0] < ann.shape[0] and pnt[1] < ann.shape[1] and pnt[2] < ann.shape[2]]
                        counts = len([cnt for cnt in counts if int(ann[tuple(cnt)]) in ids_all_structures[strt]])
                        counts_per_slice.append(counts)
                counts_per_structure.append(counts_per_slice)
            
            #now,get the total number of voxels covered by the structure in each slice
            voxels_per_structure = []
            for strt, zr in enumerate(zr_per_slice): #per structure
                voxels_per_slice = []
                for i,slc in enumerate(zr): #for each slice in the structure
                    if i < len(zr)-1: #so that the 2nd to last number in the array is considered as the beginning of the last slice
                        counts = 0 #init counter
                        iids_in_slc = np.unique(ann[slc:zr[i+1]]).astype(int)
                        iids_from_strct_in_slc = [iid for iid in iids_in_slc if iid in ids_all_structures[strt]]
                        voxels = sum(np.array([len(np.where(ann[slc:zr[i+1]] == iid)[0]) for iid in iids_from_strct_in_slc])) #total voxels of whole structure
                        voxels_per_slice.append(voxels)
                voxels_per_structure.append(voxels_per_slice)
                
            #density = count/voxelsinstructure x (scale^3)
            density_per_structure = []
            scale = 0.025 #allen
            for i, counts in enumerate(counts_per_structure):
                rho_per_slice = []
                for j,cnt in enumerate(counts):
                    rho = cnt/(voxels_per_structure[i][j]*(scale**3))
                    rho_per_slice.append(rho)
                density_per_structure.append(rho_per_slice)
            
            #make this into a dataframe now
            df = pd.DataFrame()
            df["structure"] = structures
            for i,rho_per_slice in enumerate(density_per_structure):
                for j,rho in enumerate(rho_per_slice):
                    df.loc[df.structure == structures[i], "Slice %s" % j] = rho
            
            df = df.round(2)
            df.to_csv(os.path.join(dst, "%s_75um_edge_100um_ventricle_erosion_densities.csv" % os.path.basename(src)))
            
            print ("\n Analysis Completed\n") 
        else:
            print ("\n Cells not transformed properly, check transformix (aka step 6) and try again\n")
    except:
        print("\n Path for transformed cells doesn't exist\n")
        
#%%

if __name__ == "__main__":
    
    #do the atlas-y/structure stuff outside the function
    #annotation and LUT
    ann_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_zflipped_75um_erosion_100um_ventricular_erosion.tif"
    ann_lut = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/atlas/ARA2_annotation_info_w_voxel_counts.csv"
    
    #read ann
    ann = tifffile.imread(ann_pth)
    
    #open LUT csv
    lut = pd.read_csv(ann_lut)
    
    structures = ["Caudoputamen", "corpus callosum", "Nucleus accumbens", "Hippocampal formation", "Ammon's horn", 
                  "Dentate gyrus", "Midbrain, sensory related", "Midbrain, motor related", 
                "Midbrain, behavioral state related", "Cerebellar cortex", "Cerebellar nuclei", "Primary motor area",
                "Secondary motor area", "Primary somatosensory area", "Supplemental somatosensory area", "Visual areas",
                "Cerebral cortex", "Thalamus", "Hypothalamus", "Midbrain", "Hindbrain", 
                "Cerebellum", "Striatum", "Pallidum", "fiber tracts", "Olfactory areas"]
    
    #get progeny of all large structures
    ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"
    
    with open(ontology_file) as json_file:
        ontology_dict = json.load(json_file)
    
    #find progeny for all structures
    all_structures = []
    for soi in structures:
        progeny = []
        get_progeny(dic=ontology_dict,parent_structure=soi,progeny_list=progeny)
        all_structures.append(progeny)
    
    #for each structure, find ids for themselves their progeny
    ids_all_structures = []
    for i,big_struct in enumerate(all_structures):
        id_big_struct = []
        if len(big_struct) > 0:
            for struct in big_struct:
                i_d = lut.loc[(lut.name == struct) & (lut.voxels_in_structure > 0), "id"].values
                if len(i_d) > 0: 
                    i_d = i_d[0]
                    id_big_struct.append(i_d)
        #need to add the structure itself to thte list of ids
        i_d = lut.loc[(lut.name == structures[i]) & (lut.voxels_in_structure > 0), "id"].values
        if len(i_d) > 0: 
            i_d = i_d[0]
            id_big_struct.append(i_d)
        ids_all_structures.append(id_big_struct)
    
    #find the z planes in which each of these structures + their progeny exist, get the min and max of their z plane range
    zranges = []
    for iids in ids_all_structures:
        zranges_per_iid = []
        for iid in iids:
            z,y,x = np.where(ann == iid)
            if len(z) > 0:
                zranges_per_iid.append([min(z), max(z)])
            else:
                zranges_per_iid.append([np.nan, np.nan])
        zranges.append(zranges_per_iid)
    
    #now, get the min/max z plane range per BIG STRUCTURE FROM LIST
    zr_per_struct = [(min(np.array(zr)[:, 0]), max(np.array(zr)[:, 1])) for zr in zranges]
    
    pth = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/processed"
    dst = "/jukebox/LightSheetData/pni_viral_vector_core/promoter_exp_composite_results/201912_results"    
    
    #get brainnames in fld
    brains = ['v145_none_c3_1']#os.listdir(pth)
    
    #call function
    for src in brains:
        generate_cell_density_per_slice_df(os.path.join(pth, src), dst, ann, structures,
                                           ids_all_structures, zr_per_struct)
        
    #%%
    #compond individual dfs in one
    src = os.path.join(dst, "1mo")
    csvs = [os.path.join(src, xx) for xx in os.listdir(src)]
    
    #read each, add brain column
    bigdf = []
    for i,csv in enumerate(csvs):
        df = pd.read_csv(csv, index_col = None).drop(columns = ["Unnamed: 0"])
        df["Brain"] = [os.path.basename(csv)[:-48]]*len(df) #removes added txt
        bigdf.append(df)
    
    df = pd.concat(bigdf)
    cols = ['structure', 'Brain', 'Slice 0', 'Slice 1', 'Slice 2', 'Slice 3', 'Slice 4',
       'Slice 5', 'Slice 6', 'Slice 7', 'Slice 8', 'Slice 9']
    df = df[cols]
    df = df.sort_values(by = ["structure"])
    
    df.to_csv(os.path.join(os.path.dirname(src), "1mo_composite_results_500um_slices_75um_erosion_densities.csv"))

#%%
    #plot these to make sure
    df = pd.read_csv(os.path.join(dst, "6mo_composite_results_500um_slices_75um_erosion_densities.csv")).drop(columns = ["Unnamed: 0"])
    
#    #need to plot slices as if brains - make into separate column
#    df = df.melt(id_vars = ["structure", "Brain"],
#        var_name = "Slice",
#        value_name = "Density")
#    
#    #resave them that way
#    df.to_csv(os.path.join(dst, "6mo_composite_results_500um_slices_75um_erosion_densities.csv"))
    
    #exclude the eGFP promoter
    df = df[df.Brain.str[:4] != "v190"]
    
    #categorize into lap1, lap2, ef1alpha, etc..
    df["promoter"] = ["Control"]*(len(df))
    df.loc[df.Brain.str[:3] == "v75", "promoter"] = "EF1alpha"
    df.loc[df.Brain.str[:4] == "v143", "promoter"] = "LAP1"
    df.loc[df.Brain.str[:4] == "v144", "promoter"] = "LAP2"
    df.loc[df.Brain.str[:4] == "v145", "promoter"] = "LAP1_2"
    
    #plot
    sns.barplot(x = "promoter", y = "Density", data = df[df.structure == "Hypothalamus"])
