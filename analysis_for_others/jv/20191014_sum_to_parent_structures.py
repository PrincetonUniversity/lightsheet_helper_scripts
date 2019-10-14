#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:29:36 2019

@author: wanglab
"""

import pandas as pd, os
import numpy as np

#set the id table, and annotation file paths
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"

#path to appropriate csv file
percent_density_csv_pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/cell_counts_dataframe_w_percents_density.csv"

#SET THE DESTINATION DIRECTORY HERE
dst = "/home/wanglab/Desktop"
 
#give list of structures you want to pool
pth = "/home/wanglab/Desktop/one_way_anova_all_structures.csv"

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

df = pd.read_csv(pth).drop(columns = ["Unnamed: 0"])
sois = df.name.values

#import the main excel sheet from which you will get all the voxel counts
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

#loop through and find downstream structures of a given parent
for soi in sois:
    soi = [s for s in structures if s.name==soi][0]
    print(soi.name)
    df.loc[df.name == soi.name, "voxels_in_structure"] = soi.volume_progeny
   
df.to_csv(pth, index = None)
