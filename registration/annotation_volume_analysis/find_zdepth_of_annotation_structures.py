# -*- coding: utf-8 -*-
"""
Created on Mon May 11 11:39:05 2020

@author: Zahra
"""

import os, tifffile, numpy as np, pandas as pd, json
from scipy.ndimage import center_of_mass

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
    if "msg" in list(dic.keys()): dic = dic["msg"][0]
    
    name = dic.get("name")
    children = dic.get("children")
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get("name")
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get("name")
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 
    
if __name__ == "__main__":

    #reading paths
    ann_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans.tif"
    df_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
    dst = "/jukebox/LightSheetTransfer/atlas/allen_atlas"
    ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"
    
    #get progeny of all large structures
    with open(ontology_file) as json_file:
        ontology_dict = json.load(json_file)
    
    #run if you need to get it for all structures in the lowest level of the ontology    
    ann = tifffile.imread(ann_pth)
    df = pd.read_excel(df_pth)
    
    #init z-depth column in dataframe
    df["dorsal-ventral_depth_median"] = [np.nan]*len(df)
    df["dorsal-ventral_depth_min"] = [np.nan]*len(df)
    df["dorsal-ventral_depth_max"] = [np.nan]*len(df)
    df["dorsal-ventral_depth_center_of_mass"] = [np.nan]*len(df)
    #get all structure ids
    iids = np.unique(ann).astype("float32")
    
    #iterate through ids
    for i,iid in enumerate(iids):
        #find coordinates where the id/structure exists
        zid, yid, xid = np.where(ann == iid)
        #take the median of the z-coordinate to get 'center' in z
        x_median = np.median(xid)
        x_min = np.min(xid)
        x_max = np.max(xid)
        #find the center of mass
        ann_iid = np.zeros_like(ann)
        ann_iid[zid,yid,xid] = 1
        zcm, ycm, xcm = center_of_mass(ann_iid)
        if i%50==0:
            print("******%s******\n" % iid)
            print("******median dorsal-ventral_coordinate is %s******\n" % x_median)
            print("******center of mass is %s******\n" % xcm)
        df.loc[df.id == iid, "dorsal-ventral_depth_median"] = x_median
        df.loc[df.id == iid, "dorsal-ventral_depth_min"] = x_min
        df.loc[df.id == iid, "dorsal-ventral_depth_max"] = x_max
        df.loc[df.id == iid, "dorsal-ventral_depth_center_of_mass"] = np.round(xcm, 1)
    
    print("******saving to dataframe...******\n")
    df.to_excel(os.path.join(dst, "allen_id_table_w_voxel_counts_n_DVcoord.xlsx"))
    
    #%%
    
    # df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_n_DVcoord.csv"
    # df = pd.read_csv(df_pth)
    
    # str_pth = "/jukebox/wang/Jess/lightsheet_output/structuresorting/structures_4.13.20.csv"
    # structs_df = pd.read_csv(str_pth, index_col = None, header = None)
    # structs = [xx[0] for xx in np.array(structs_df)]
    
    # #init column to put z-depth
    # structs_df["dorsal-ventral_coordinate"] = [np.nan]*len(structs_df)
    # #iterate through structs and find median z-depth of their child structures
    # for struct in structs:
    #     print("******%s******\n" % struct)
    #     progeny = []; z_depths = []
    #     get_progeny(ontology_dict, struct, progeny)
    #     #add original structure z-depth to list
    #     if struct in df.name.values: z_depths.append(df.loc[df.name == struct, 
    #                         "dorsal-ventral_coordinate"].values[0])
    #     for progen in progeny:
    #         if progen in df.name.values:
    #             z_depths.append(df.loc[df.name == progen, "dorsal-ventral_coordinate"].values[0])
    #     structs_df.loc[structs_df[0] == struct, "dorsal-ventral_coordinate"] = np.nanmedian(np.array(z_depths))
        
    # #format and export
    # structs_df.columns = ["name", "dorsal-ventral_coordinate"]
    # structs_df.to_excel(os.path.join(os.path.dirname(str_pth), "structures_4.13.20_w_DVcoord.xlsx"), index = None)
