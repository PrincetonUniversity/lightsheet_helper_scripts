#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:24:45 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl, json, itertools, statsmodels.api as sm
from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.registration.register import change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files
from tools.registration.transform_list_of_points import point_transformix, unpack_pnts
from tools.utils.io import makedir
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass
import matplotlib.colors
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop
 
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"
fig_dst = "/home/wanglab/Desktop"

ann_pth = os.path.join(dst, "atlases/sagittal_allen_ann_25um_iso_60um_edge_80um_ventricular_erosion.tif")
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"

#cut annotation file in middle
ann = tifffile.imread(ann_pth)
plt.imshow(ann[300])
z,y,x = ann.shape
#make sure each halves are same dimension as original ann
ann_left = np.zeros_like(ann)
ann_left[:int(z/2), :, :] = ann[:int(z/2), :, :] #cut in the middle in x
ann_right = np.zeros_like(ann)
ann_right[int(z/2):, :, :] = ann[int(z/2):, :, :]
plt.imshow(ann_left[120])

#collect 
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20170410_tp_bl6_lob6a_ml_repro_01", "20160823_tp_bl6_cri_500r_02", "20180417_jg59_bl6_cri_03",
"20170207_db_bl6_crii_1300r_02", "20160622_db_bl6_unk_01", "20161205_tp_bl6_sim_750r_03",
"20180410_jg51_bl6_lob6b_04", "20170419_db_bl6_cri_rpv_53hr", "20170116_tp_bl6_lob6b_lpv_07",
"20170411_db_bl6_crii_mid_53hr", "20160822_tp_bl6_crii_1500r_06", "20160920_tp_bl6_lob7_500r_03",
"20170207_db_bl6_crii_rpv_01", "20161205_tp_bl6_sim_250r_02", "20161207_db_bl6_lob6a_500r_53hr",
"20170130_tp_bl6_sim_rlat_05", "20170115_tp_bl6_lob6b_500r_05", "20170419_db_bl6_cri_mid_53hr",
"20161207_db_bl6_lob6a_850r_53hr", "20160622_db_bl6_crii_52hr_01", "20161207_db_bl6_lob6a_50rml_53d5hr",
"20161205_tp_bl6_lob45_1000r_01", "20160801_db_l7_cri_01_mid_64hr"]    

src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites"

imgs = [os.path.join(src, xx+".tif.tif") for xx in brains]

#pool brain names and L/R designation into dict
lr_dist = {}
inj_vox = {}

#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    print(brain)
    inj_vol = tifffile.imread(img)
    z,y,x = inj_vol.shape
    
    z_c,y_c,x_c = center_of_mass(inj_vol)
    #take distance from center to arbitrary "midline" (aka half of z axis)
    dist = z_c-(z/2)
    #save to dict 
    lr_dist[brain[:-8]] = dist
    inj_vox[brain[:-8]] = inj_vol
    
    if dist < 0:
        print("brain {} has a left-sided injection\n".format(brain))
    elif dist > 0:
        print("brain {} has a right-sided injection\n".format(brain))
    else:
        print("brain has an injection close to midline so not considering it rn\n")


#make structures
#get injection fractions
inj_raw = np.array([inj.astype(bool) for inj in inj_vox.values()])
    
atl_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif")[:, 423:, :] #cropping coz tom segmented this
ann_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")[:, 423:, :]
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)
     
#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lingula (I)": 912,
        "Lobule II": 976,
        "Lobule III": 984,
        "Lobule IV-V": 1091,
        "Lobule VIa": 936,
        "Lobule VIb": 1134,
        "Lobule VII": 944,
        "Lobule VIII": 951,
        "Lobule IX": 957, #uvula IX
        "Lobule X": 968, #nodulus X
        "Simplex lobule": 1007, #simplex
        "Crus 1": 1056, #crus 1
        "Crus 2": 1064, #crus 2
        "Paramedian lobule": 1025, #paramedian lob
        "Copula pyramidis": 1033 #copula pyramidis
        }
ak = np.asarray([k for k,v in iids.items()])

atlas_rois = {}
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    ann_blank = np.zeros_like(ann_raw)
    ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
    atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary

expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

#pooled injections
ak_pool = np.array(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X", #no simpplex injections
                 "Simplex", "Crus I", "Crus II", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
#get n's after pooling
primary_lob_n = np.array([len(np.where(primary_pool == i)[0]) for i in range(max(primary_pool)+1)])
#normalization  of inj site
frac_of_inj_pool_norm = np.asarray([brain/brain.sum() for brain in frac_of_inj_pool])

lr_brains = list(lr_dist.keys())
atl_dst = os.path.join(dst, "pma_to_aba"); makedir(atl_dst)
id_table = pd.read_excel(df_pth)
#%%
#------------------------------------------------------------------------------------------------------------------------------
#NOTE THAT ONLY HAVE TO DO THIS ONCE!!!! DO NOT NEED TO DO AGAIN UNLESS DOUBLE CHECKIHG
#transform points to allen atlas space
nc_lst = brains

#get brains that we actually need to get cell counts from
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_thalamus/transformed_points"
post_transformed = [os.path.join(src, os.path.join(xx, "transformed_points/posttransformed_zyx_voxels.npy")) for xx in lr_brains]
transformfiles = ["/jukebox/wang/zahra/aba_to_pma/TransformParameters.0.txt",
                  "/jukebox/wang/zahra/aba_to_pma/TransformParameters.1.txt"]

##collect 
#for fl in post_transformed:
#    arr = np.load(fl)
#    #make into transformix-friendly text file
#    brain = os.path.basename(os.path.dirname(os.path.dirname(fl)))
#    print(brain)
#    transformed_dst = os.path.join(atl_dst, brain); makedir(atl_dst)
#    pretransform_text_file = create_text_file_for_elastix(arr, transformed_dst)
#        
#    #copy over elastix files
#    trfm_fl = modify_transform_files(transformfiles, transformed_dst) 
#    change_transform_parameter_initial_transform(trfm_fl[0], 'NoInitialTransform')
#   
#    #run transformix on points
#    points_file = point_transformix(pretransform_text_file, trfm_fl[-1], transformed_dst)
#    
#    #convert registered points into structure counts
#    converted_points = unpack_pnts(points_file, transformed_dst) 
#    
##%%
##------------------------------------------------------------------------------------------------------------------------------    
#def transformed_cells_to_allen(fld, ann, dst, fl_nm):
#    """ consolidating to one function bc then no need to copy/paste """
#    dct = {}
#    
#    for fl in fld:
#        converted_points = os.path.join(fl, "posttransformed_zyx_voxels.npy")
#        print(converted_points)
#        point_lst = transformed_pnts_to_allen_helper_func(np.load(converted_points), ann, order = "ZYX")
#        df = count_structure_lister(id_table, *point_lst).fillna(0)
#        #for some reason duplicating columns, so use this
#        nm_cnt = pd.Series(df.cell_count.values, df.name.values).to_dict()
#        fl_name = os.path.basename(fl)
#        dct[fl_name]= nm_cnt
#        
#    #unpack
#    index = dct[list(dct.keys())[0]].keys()
#    columns = dct.keys()
#    data = np.asarray([[dct[col][idx] for idx in index] for col in columns])
#    df = pd.DataFrame(data.T, columns=columns, index=index)
#    
#    #save before adding projeny counts at each level
#    df.to_pickle(os.path.join(dst, fl_nm))
#    
#    return os.path.join(dst, fl_nm)
#
#pma2aba_transformed = [os.path.join(atl_dst, xx) for xx in lr_brains]
##collect counts from right side
#right = transformed_cells_to_allen(pma2aba_transformed, ann_right, dst, "thal_right_side_no_prog_at_each_level_allen_atl.p")
##collect counts from left side
#left = transformed_cells_to_allen(pma2aba_transformed, ann_left, dst, "thal_left_side_no_prog_at_each_level_allen_atl.p")

#import dict of cells by region
r_cells_regions = pckl.load(open(os.path.join(dst, "thal_right_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
r_cells_regions = r_cells_regions.to_dict(orient = "dict")      

contra = {}; ipsi = {} #collect contra and ipsi frame
for k,v in r_cells_regions.items():
    if lr_dist[k] < 0:
        contra[k] = v
    else:
        ipsi[k] = v

#LEFT SIDE
l_cells_regions = pckl.load(open(os.path.join(dst, "thal_left_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
l_cells_regions = l_cells_regions.to_dict(orient = "dict")      

for k,v in l_cells_regions.items():
    if lr_dist[k] > 0:
        contra[k] = v
    else:
        ipsi[k] = v
        
contra_df = pd.DataFrame(contra)
contra_df.to_csv(os.path.join(dst, "data/thal_contra_counts_23_brains.csv")) 

ipsi_df = pd.DataFrame(ipsi)
ipsi_df.to_csv(os.path.join(dst, "data/thal_ipsi_counts_23_brains.csv"))         

#%%
cells_regions_pth = os.path.join(dst, "data/thal_contra_counts_23_brains.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.025
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

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

#get progeny of all large structures
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for all of neocortex
sois = ["Thalamus", "Ventral posteromedial nucleus of the thalamus", "Ventral posterolateral nucleus of the thalamus",
          "Ventral anterior-lateral complex of the thalamus", "Ventral medial nucleus of the thalamus", "Anteroventral nucleus of thalamus", 
          "Reticular nucleus of the thalamus", "Ventral part of the lateral geniculate complex", "Mediodorsal nucleus of thalamus",
          "Submedial nucleus of the thalamus", "Nucleus of reuniens", "Paraventricular nucleus of the thalamus", 
          "Central lateral nucleus of the thalamus", "Parafascicular nucleus", "Posterior complex of the thalamus",
          "Lateral dorsal nucleus of thalamus", "Lateral posterior nucleus of the thalamus", "Lateral habenula"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)



#%%
