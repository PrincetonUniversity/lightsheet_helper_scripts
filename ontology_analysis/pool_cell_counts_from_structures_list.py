#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 14:26:37 2019

@author: wanglab
"""

import pandas as pd, os
import numpy as np

#set the id table, and annotation file paths
df_pth = "/home/wanglab/mounts/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ann_pth = "/home/wanglab/mounts/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_60um_erosion.tif"
ontology_file = "/path/to/.json file"
#path to appropriate csv file
percent_density_csv_pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/cell_counts_dataframe_w_percents_density.csv"

#SET THE DESTINATION DIRECTORY HERE
dst = "/home/wanglab/Desktop"
 
#give list of structures you want to pool
pth = "/home/wanglab/mounts/wang/Jess/lightsheet_output/201904_ymaze_cfos/structures.csv"

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
    name = dic.get('name')
    children = dic.get('children')
    if name == input_nodename:
        for child in children: # child is a dict
            child_name = child.get('name')
            progeny_list.append(child_name)
            get_progeny(child,input_nodename=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get('name')
        get_progeny(child,input_nodename=input_nodename,progeny_list=progeny_list)
    return 

if __name__ == '__main__':

    """ pools regions together based on allen name """    

    """ first load in the json dictionary containing the ontology """
    with open(ontology_file) as json_file:
        ontology_dict = json.load(json_file)

    sois = pd.read_csv(pth)
    sois = [xx[0] for xx in sois.values]
       
    #set variables
    orgdf = pd.read_csv(percent_density_csv_pth)

    #init
    tdf = orgdf.copy()

    #we have to do this per animal
    #get the brains
    brains = tdf.Brain.unique()

    #make a big dicitionary to put it all together
    pooled_counts_an = {}

    #iterate through the conditions
    for an in brains[0:1]:
        df = tdf[tdf.Brain == an]
        condition = df["Condition"].unique()
        pooled_counts = {}
        
        #loop through and find downstream structures of a given parent
        for soi in sois[0:1]:
            # soi = [s for s in structures if s.name==soi][0]
            # print(soi.name)
            # progeny = [str(xx.name) for xx in soi.progeny]
            progeny = []
            progeny = get_progeny(dic=ontology_dict,parent_structure=soi,progeny_list=progeny)
            counts = [] #store counts in this list
            val = df.loc[(df.name == soi.name), "percent"].values
            if val.shape[0] > 0: counts.append(val[0])
            if len(progeny) > 0:
                for progen in progeny:
                    val = df.loc[(df.name == progen), "percent"].values
                    if val.shape[0] > 0: counts.append(val[0])
                    
            #sum counts in parent structures            
            pooled_counts[soi.name] = np.sum(np.asarray(counts))
        
        #fill other details and add to big dict
        pooled_counts_an[an] = pooled_counts
        pooled_counts["group"] = condition[0]
        
    #make into giant dataframe
    main_df = pd.DataFrame.from_dict(pooled_counts_an, orient = "index")
    #sort by group so makes more sense
    main_df = main_df.sort_values(by=["group"])
    main_df.index.name = "animal"

    main_df.to_csv(os.path.join(dst, "select_structures_percent_counts_for_visualization.csv"))

    import itertools

    rotate_df = pd.DataFrame()
    struct = [list(itertools.repeat(xx, 33)) for xx in main_df.columns.values[:-1]]
    struct = pd.Series(list(itertools.chain.from_iterable(struct)))
    rotate_df["name"] = struct

    vals = [pd.Series(main_df[xx].values) for xx in main_df.columns.values[:-1]]    
    rotate_df["percent"] = pd.concat(vals, ignore_index = True)

    ans = list(itertools.repeat(main_df.index.values, len(struct)))
    ans = pd.Series(list(itertools.chain.from_iterable(ans)))
    rotate_df["animal"] = ans

    groups = list(itertools.repeat(main_df.group.values, len(struct)))
    groups = pd.Series(list(itertools.chain.from_iterable(groups)))
    rotate_df["condition"] = groups

    #save
    rotate_df.to_csv(os.path.join(dst, "select_structures_percent_counts_for_plots.csv"), index = False)

    print("saved in :{}".format(dst))
    

