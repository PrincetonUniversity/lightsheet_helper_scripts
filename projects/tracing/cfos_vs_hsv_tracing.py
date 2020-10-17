#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:45:34 2020

@author: wanglab
"""
import os, numpy as np, pandas as pd, scipy, itertools, sys,json

src="/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2020_mapping_paper"
brains = ["/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk01",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk02",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk03",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk04",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk05",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk06",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk07",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk08",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk10",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk11",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp02",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp06",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp08",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp09",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files_reim/201701_tp01",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files_reim/201701_tp07",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files_reim/201701_tpal",
 "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files_reim/201701_tpbe"]

#get files
lst = [os.path.join(br,"Annotated_counts_Allen_20201016.csv") for br in brains]
#conditions
nms = [os.path.basename(br) for br in brains]

cond = ["stim","stim","stim","stim","stim","stim","stim","stim","stim","stim","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl"]
conditions = {n:c for n,c in zip(nms, cond)}

#set appropriate paths
pth = os.path.join(src, "cfos_cell_counts_Allen_20201016.csv")
df_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
ann_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans.tif"
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

def generate_data_frame(conditions, lst, pth, flnm):
    """ 
    used to make a pooled csv file of all cell counts in an experiment
    inputs:
        conditions: zip of file names + condition
        lst: list of file names run through analysis
        pth: path to save csv output
    """
    #generate data frame
    bglst=[]
    for fl in lst:
        #extract out info
        nm = os.path.basename(os.path.dirname(fl))
        #make dataframe
        df = pd.read_csv(fl)[1:] #remove previous headers
        print(nm, df.shape)
        df = df.replace(np.nan, "", regex=True)
        df["Brain"] = nm
        df["Condition"] = conditions[nm]
        bglst.append(df)
    
    df = pd.concat(bglst)
    df["counts"] = df["counts"].apply(int)

    df.drop(columns = ["Unnamed: 0"]).to_csv(pth, index = None)
    
    return pth

#run
csv_pth = generate_data_frame(conditions, lst, pth, "Annotated_counts_Allen_20201016")

def generate_percent_counts_and_density_per_region(src, csv_pth):
    """ generates another column in the dataframe that is just # counts / total counts in brain and density per region"""    
    #read csv
    df = pd.read_csv(csv_pth)
    
    #set resolution of atlas
    scale = 0.025 ##mm/voxel
    
    #get all brain names
    brains = np.unique(df.Brain.values)
    
    #save total counts in dict
    total_counts = {}
    percents = {}
    #for each brain, get total counts
    for brain in brains:
        total_counts[brain] = df[df.Brain == brain].counts.sum(0)    
           
    percents = [df[df.Brain == brain].counts.apply(lambda x: (x/total_counts[brain])*100).astype("float64") for brain in brains]
    
    #concantenate together
    df_percent = pd.concat(percents)
    
    df["percent"] = df_percent
    df["density"] = df[df.voxels_in_structure > 0].apply(lambda x:x.counts/(float(x.voxels_in_structure*(scale**3))), 1)
    
    #export
    df.to_csv(pth,index=None)
    
    return os.path.join(pth)

#run
percent_density_csv_pth = generate_percent_counts_and_density_per_region(src, csv_pth)

#%%
#read
cells_regions=pd.read_csv(percent_density_csv_pth,index_col=None)

#pool nc structures and generate dataframe
def get_progeny(dic,parent_structure,progeny_list):
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

#get progeny of all large structures
with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for all of neocortex
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        br=[] #per brain
        try:
            for nm in nms:
                try:
                    br.append(cells_regions[cells_regions.name==progen].loc[cells_regions.Brain == nm, "counts"].values[0])
                except:
                    br.append(0)
            counts.append(br)
        except:
            counts.append([0]*len(nms))
    counts_per_struct.append(np.array(counts).sum(axis = 0))
    
counts_per_struct = np.array(counts_per_struct)

ann_df=pd.read_excel(df_pth)
#get volumes
vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        
scale_factor=0.025 #allen
density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T

#layer p counts maps
pcounts = np.array([xx/sum(xx) for xx in counts_per_struct.T])*100

sois = np.array(sois)[np.argsort(vol)]
pcounts = pcounts.T[np.argsort(vol)].T
density = density.T[np.argsort(vol)].T
counts_per_struct = counts_per_struct[np.argsort(vol)].T

#save out density
df=pd.DataFrame(density)
df.index=nms
df.columns=sois
df["condition"]=cond
df=df.round(2)
df.to_csv(os.path.join(src,"cfos_density_NC_only_pooled_Allen_20201016.csv"))

#save out counts
df=pd.DataFrame(counts_per_struct)
df.index=nms
df.columns=sois
df["condition"]=cond
df=df.round(2)
df.to_csv(os.path.join(src,"cfos_cell_counts_NC_only_pooled_Allen_20201016.csv"))

#transform cells to PMA and re-run analysis
from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.registration.register import change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files
from tools.registration.transform_list_of_points import point_transformix, unpack_pnts

transformfiles = ["/jukebox/wang/zahra/aba_to_pma/TransformParameters.0.txt",
                  "/jukebox/wang/zahra/aba_to_pma/TransformParameters.1.txt"]
atl_dst="/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/aba_to_pma_transformed_cells"
#collect 
for br in brains:
    arr = np.load(os.path.join(br,"clearmap_cluster_output/cells_transformed_to_Atlas.npy"))
    #change to zyx!!!
    arr=np.array([arr[:,2],arr[:,1],arr[:,0]]).T
    #make into transformix-friendly text file
    transformed_dst = os.path.join(atl_dst, os.path.basename(br)); 
    if not os.path.exists(transformed_dst): os.mkdir(transformed_dst)
    pretransform_text_file = create_text_file_for_elastix(arr, transformed_dst)
    #copy over elastix files
    trfm_fl = modify_transform_files(transformfiles, transformed_dst) 
    change_transform_parameter_initial_transform(trfm_fl[0], 'NoInitialTransform')
    #run transformix on points
    points_file = point_transformix(pretransform_text_file, trfm_fl[-1], transformed_dst)
    #convert registered points into structure counts
    converted_points = unpack_pnts(points_file, transformed_dst) 

#make into dataframe
def transformed_cells_to_allen(fld, ann, dst, fl_nm, id_table=df_pth):
    """ consolidating to one function bc then no need to copy/paste """
    dct = {}
    for fl in fld:
        converted_points = os.path.join(fl, "posttransformed_zyx_voxels.npy")
        print(converted_points)
        point_lst = transformed_pnts_to_allen_helper_func(np.load(converted_points), ann, order="ZYX")
        df = count_structure_lister(id_table, *point_lst).fillna(0)
        #for some reason duplicating columns, so use this
        nm_cnt = pd.Series(df.cell_count.values, df.name.values).to_dict()
        fl_name = os.path.basename(fl)
        dct[fl_name]= nm_cnt
    #unpack
    index = dct[list(dct.keys())[0]].keys()
    columns = dct.keys()
    data = np.asarray([[dct[col][idx] for idx in index] for col in columns])
    df = pd.DataFrame(data.T, columns=columns, index=index)
    #save before adding projeny counts at each level
    df.to_csv(os.path.join(dst, fl_nm))
    return os.path.join(dst, fl_nm)
#run
pma2aba_transformed = [os.path.join(atl_dst, os.path.basename(xx)) for xx in brains]
df = transformed_cells_to_allen(pma2aba_transformed, tifffile.imread(ann_pth), src, "cfos_cell_counts_NC_only_pooled_PMA_20201016.csv")

#recalculate density for PMA cells
#first calculate counts across entire nc region
cells_regions = pd.read_csv(df)
cells_regions["name"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns=["Unnamed: 0"])
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    counts.append([cells_regions.loc[cells_regions.name == soi, nm].values[0] for nm in nms])
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.name == progen, nm].values[0] for nm in nms])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)
    
#change to the PMA
df_pth="/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ann_df=pd.read_excel(df_pth)
#get volumes
vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        
scale_factor=0.020 #pma
density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T

#layer p counts maps
pcounts = np.array([xx/sum(xx) for xx in counts_per_struct.T])*100

sois = np.array(sois)[np.argsort(vol)]
pcounts = pcounts.T[np.argsort(vol)].T
density = density.T[np.argsort(vol)].T

#save out density
df=pd.DataFrame(density)
df.index=nms
df.columns=sois
df["condition"]=cond
df=df.round(2)
df.to_csv(os.path.join(src,"cfos_density_NC_only_pooled_PMA_20201016.csv"))

#mapping to check
cells=np.load("/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/aba_to_pma_transformed_cells/201701_mk08/posttransformed_zyx_voxels.npy")
ann=tifffile.imread(ann_pth)
cellmap=np.zeros(ann.shape)
for cell in cells.astype(int):
    cellmap[cell[0],cell[1],cell[2]]=1
tifffile.imsave("/home/wanglab/Desktop/test.tif",cellmap.astype("uint8"))
