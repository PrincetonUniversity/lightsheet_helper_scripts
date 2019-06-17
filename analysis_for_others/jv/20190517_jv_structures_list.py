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

#path to appropriate csv file
percent_density_csv_pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/cell_counts_dataframe_w_percents_density.csv"

#SET THE DESTINATION DIRECTORY HERE
dst = "/home/wanglab/Desktop"
 
#give list of structures you want to pool
pth = "/home/wanglab/mounts/wang/Jess/lightsheet_output/201904_ymaze_cfos/structures.csv"

#run this first cell the way it is, imports tom"s class for allen ontology
class structure:
    """Class to represent a brain structure
    """
    kind="structure"    
    def __init__(self, idnum, excelfl, units="pixels", scale_factor=1):
        self.idnum=float(idnum)        #id of structure 
        self.idnum_int = int(idnum)    #id of structure as float, sometimes useful for comparisions
        self.excelfl=excelfl            #path to excelfl
        self.acronym=()                 #acronym of structure        
        self.name=()                    #allen name of strcutre
        self.cellcount=()               #determined count of actual structure
        self.cellcount_progeny=()       #determined count of progeny + structure
        self.parent=()                  #parent id, acronym, name
        self.children=[]                #children one level down (first sublevel)
        self.progeny=[]                 #all children, grandchildren, etc
        self.progeny_pixels=[idnum]     #all progeny pixel IDs and it"s own
        self.progenitor_chain=[]        #list moving up the heirarchy. i.e.: LGN->Thal->Interbrain
        self.volume=()                  #number of voxels, scale factor is not accounted for unless provided
        self.volume_progeny=()           #number of voxels, scale factor is not accounted for unless provided of progeny+structure
        self.units=units                #units that volume is displayed as
        self.scale_factor = scale_factor#unittopixel conversion. I.e. aba=25um/pixels == 25
        
    def add_name(self, nm): 
        self.name=nm
    def add_acronym(self, acronym): 
        self.acronym=acronym
    def add_cellcount(self, num): 
        self.cellcount=num
    def add_cellcount_progeny(self, num): 
        self.cellcount_progeny=num
    def add_parent(self, nm): 
        self.parent=nm
    def add_child(self, nm): 
        self.children.append(nm)
    def add_progeny(self, nm): 
        self.progeny.append(nm)
    def add_progeny_pixels(self, nm): 
        self.progeny_pixels.append(nm)
    def create_progenitor_chain(self, lst):
        self.progenitor_chain=lst

def find_progeny(struct, df):
    #find children (first sublevel)   
    #children = [x for x in df[df["parent_structure_id"] == str(struct.idnum)].itertuples()]
    children = [x for x in df[df["parent_structure_id"] == struct.idnum].itertuples()]
    
    #find all progeny
    allchildstructures=[]
    while len(children) > 0:
        child = children.pop()
        allchildstructures.append(child)
        #kiditems = df[df["parent_structure_id"] == str(child.id)]
        kiditems = df[df["parent_structure_id"] == child.id]
        for kid in kiditems.itertuples():
            allchildstructures.append(kid)            
            #if kid has children append them list to walk through it
            #if len(df[df["parent_structure_id"] == str(kid.id)]) != 0:
            if len(df[df["parent_structure_id"] == kid.id]) != 0:                
                children.append(kid)
    #remove duplicates
    allchildstructures = list(set(allchildstructures))
    #add progeny to structure_class
    [struct.add_progeny(xx) for xx in allchildstructures]
    #add_progeny_pixels to structure_class
    [struct.add_progeny_pixels(xx.id) for xx in allchildstructures] #<---11/1/17, this was atlas_id, changing to id
    #add progeny count
    struct.add_cellcount_progeny(struct.cellcount + sum([int(xx.cell_count) for xx in allchildstructures]))
    if "voxels_in_structure" in df.columns: struct.volume_progeny = struct.volume + sum([int(xx.voxels_in_structure) for xx in allchildstructures])
    
    return struct
    
def create_progenitor_chain(structures, df, verbose = False):

    new_structures=[]
    for struct in structures:
        if verbose: print(struct.name)
        if struct.name == "root" or struct.name == "Basic cell groups and regions":
            pass
        else:
            chain = []
            loop = True
            current_struct = struct
            while loop:
                #find parent
                parent = current_struct.parent[1]
                if parent == "nan" or parent == "null" or parent == "root": break
                else:#append
                    chain.append(parent)
                    current_struct = [xx for xx in structures if xx.name == parent][0]
            struct.create_progenitor_chain(chain)
            new_structures.append(struct)    

    return new_structures

def make_structure_objects(excelfl, remove_childless_structures_not_repsented_in_ABA = False, ann_pth = None, verbose = False):

    #load
    df = pd.read_excel(excelfl)
    if not "cell_count" in df.columns: df["cell_count"] = 0
    #make structure objects for each df row
    structures = []
    
    for row in df.itertuples():
        #basic structure class features
        struct = structure(row.id , excelfl) #make object using ID and location of excelfl
        struct.add_name(str(row.name)) #add name        
        struct.add_acronym(str(row.acronym)) #add acronym
        struct.add_cellcount(row.cell_count) #add cell count
        struct.add_parent((row.parent_structure_id, str(row.parent_name), str(row.parent_acronym))) #parent id, acronym, name
        if "voxels_in_structure" in df.columns: struct.volume = row.voxels_in_structure
        
        #find children (first sublevel)   
        #children = [x for x in df[df["parent_structure_id"] == str(struct.idnum_int)].itertuples()]
        children = [x for x in df[df["parent_structure_id"] == struct.idnum_int].itertuples()]
        struct.add_child([xx for xx in children])
              
        #add structure to structures list
        structures.append(struct)
    #find progeny (all sublevels)
    structures = [find_progeny(struct, df) for struct in structures]
    
    #create progenitor_chain
    structures = create_progenitor_chain(structures, df)
    ###        
    return structures
 
#build structures class
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#%%

""" pools regions together based on allen name """    

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
for an in brains:
    df = tdf[tdf.Brain == an]
    condition = df["Condition"].unique()
    pooled_counts = {}
    
    #loop through and find downstream structures of a given parent
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print(soi.name)
        progeny = [str(xx.name) for xx in soi.progeny]
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

#rotate the df to make it eaiser to plot things
rotate_df = pd.DataFrame()
structures = pd.Series(main_df.columns.values[1:])
rotate_df["name"] = pd.concat([structures]*33)
vals = [pd.Series(main_df[xx].values) for xx in main_df.columns.values[1:]]    
rotate_df["percent"] = pd.concat(vals, ignore_index = True)
ans = [pd.Series([xx]*structures.shape[0]) for xx in main_df.index.values]
rotate_df["animal"] = pd.concat(ans)
groups = [pd.Series([xx]*structures.shape[0]) for xx in main_df.group.values]
rotate_df["condition"] = pd.concat(groups)

#save
rotate_df.to_csv(os.path.join(dst, "select_structures_percent_counts_for_plots.csv"), index = False)

print("saved in :{}".format(dst))
    

