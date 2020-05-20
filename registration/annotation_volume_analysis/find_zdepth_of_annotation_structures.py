# -*- coding: utf-8 -*-
"""
Created on Mon May 11 11:39:05 2020

@author: Zahra
"""

import os, tifffile, numpy as np, pandas as pd, json

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
    ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
    df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
    dst = "/jukebox/LightSheetTransfer/atlas"
    ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"
    
    #get progeny of all large structures
    with open(ontology_file) as json_file:
        ontology_dict = json.load(json_file)
    
    #run if you need to get it for all structures in the lowest level of the ontology    
    # ann = tifffile.imread(ann_pth)
    # df = pd.read_excel(df_pth)
    
    # #init z-depth column in dataframe
    # df["dorsal-ventral_coordinate"] = [np.nan]*len(df)
    # #get all structure ids
    # iids = np.unique(ann).astype("float32")
    
    # #iterate through ids
    # for i,iid in enumerate(iids):
    #     #find coordinates where the id/structure exists
    #     zid, yid, xid = np.where(ann == iid)
    #     #take the median of the z-coordinate to get 'center' in z
    #     x_median = np.median(xid)
    #     if i%50==0:
    #         print("******%s******\n" % iid)
    #         print("******median dorsal-ventral_coordinate is %s******\n" % x_median)
    #     df.loc[df.id == iid, "dorsal-ventral_coordinate"] = x_median
    
    # print("******saving to dataframe...******\n")
    # df.to_csv(os.path.join(dst, "ls_id_table_w_voxelcounts_n_DVcoord.csv"))
    
    #%%
    
    df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_n_DVcoord.csv"
    df = pd.read_csv(df_pth)
    
    str_pth = "/jukebox/wang/Jess/lightsheet_output/structuresorting/structures_4.13.20.csv"
    structs_df = pd.read_csv(str_pth, index_col = None, header = None)
    structs = [xx[0] for xx in np.array(structs_df)]
    
    #init column to put z-depth
    structs_df["dorsal-ventral_coordinate"] = [np.nan]*len(structs_df)
    #iterate through structs and find median z-depth of their child structures
    for struct in structs:
        print("******%s******\n" % struct)
        progeny = []; z_depths = []
        get_progeny(ontology_dict, struct, progeny)
        #add original structure z-depth to list
        if struct in df.name.values: z_depths.append(df.loc[df.name == struct, 
                            "dorsal-ventral_coordinate"].values[0])
        for progen in progeny:
            if progen in df.name.values:
                z_depths.append(df.loc[df.name == progen, "dorsal-ventral_coordinate"].values[0])
        structs_df.loc[structs_df[0] == struct, "dorsal-ventral_coordinate"] = np.nanmedian(np.array(z_depths))
        
    #format and export
    structs_df.columns = ["name", "dorsal-ventral_coordinate"]
    structs_df.to_excel(os.path.join(os.path.dirname(str_pth), "structures_4.13.20_w_DVcoord.xlsx"), index = None)
