#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 15:44:38 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl, statsmodels.api as sm, itertools
from skimage.external import tifffile
import matplotlib.colors as colors, json
from scipy.ndimage.measurements import center_of_mass

from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister, change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files, point_transformix, unpack_pnts
from tools.utils.io import makedir
from tools.analysis.network_analysis import make_structure_objects

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#set paths, variables, etc.
dst = "/jukebox/wang/zahra/tracing_projects/prv/"
fig_dst = "/home/wanglab/Desktop"
pma_ann_pth = os.path.join(dst, "pma_annotation_sagittal_atlas_20um_iso_60um_edge_160um_vntric_erosion.tif")
#cut annotation file in middle
ann = tifffile.imread(pma_ann_pth)
plt.imshow(ann[300])
z,y,x = ann.shape
#make sure each halves are same dimension as original ann
ann_left = np.zeros_like(ann)
ann_left[:int(z/2), :, :] = ann[:int(z/2), :, :] #cut in the middle in z
ann_right = np.zeros_like(ann)
ann_right[int(z/2):, :, :] = ann[int(z/2):, :, :]
plt.imshow(ann_left[120])

#collect 
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20180205_jg_bl6f_prv_01", "20180205_jg_bl6f_prv_02", "20180205_jg_bl6f_prv_03", "20180205_jg_bl6f_prv_04", 
          "20180215_jg_bl6f_prv_05", "20180215_jg_bl6f_prv_06", "20180215_jg_bl6f_prv_08", "20180215_jg_bl6f_prv_09", 
           "20180305_jg_bl6f_prv_11", "20180305_jg_bl6f_prv_12", "20180305_jg_bl6f_prv_13","20180306_jg_bl6f_prv_14", 
           "20180305_jg_bl6f_prv_15", "20180312_jg_bl6f_prv_17", "20180326_jg_bl6f_prv_37",
           "20180313_jg_bl6f_prv_21", "20180313_jg_bl6f_prv_23", "20180313_jg_bl6f_prv_24", "20180313_jg_bl6f_prv_25",
           "20180322_jg_bl6f_prv_27", "20180322_jg_bl6f_prv_28", "20180323_jg_bl6f_prv_30", "20180326_jg_bl6f_prv_33", 
           "20180326_jg_bl6f_prv_34", "20180326_jg_bl6f_prv_35"]
    
inj_src = os.path.join(dst, "prv_injection_sites")
imgs = [os.path.join(inj_src, xx+".tif") for xx in brains]

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
    lr_dist[brain[:-4]] = dist
    inj_vox[brain[:-4]] = inj_vol
    
    if dist < 0:
        print("brain {} has a left-sided injection\n".format(brain))
    elif dist > 0:
        print("brain {} has a right-sided injection\n".format(brain))
    else:
        print("brain has an injection close to midline so not considering it rn\n")


#get injection fractions
inj_raw = np.array([inj.astype(bool) for inj in inj_vox.values()])
    
atl_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif")
ann_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")
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
                 "Crus I", "Crus II", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
#get n's after pooling
primary_lob_n = np.array([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])

#make structures
#FIXME: for some reason the allen table does not work on this, is it ok to use PMA        
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=pma_ann_pth)

#%%
#set variables
lr_brains = list(lr_dist.keys())

trnsfrm_dst = os.path.join(dst, "prv_transformed_cells")
id_table = pd.read_excel(df_pth)

#post_transformed = [os.path.join(trnsfrm_dst, os.path.join(xx, "transformed_points/posttransformed_zyx_voxels.npy")) for xx in lr_brains]
#transformfiles = ["/jukebox/wang/zahra/aba_to_pma/TransformParameters.0.txt",
#                  "/jukebox/wang/zahra/aba_to_pma/TransformParameters.1.txt"]
#########################################NO NEED TO RUN AGAIN IF ALREADY RUN ONCE################################################################
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
#    change_transform_parameter_initial_transform(trfm_fl[0], "NoInitialTransform")
#   
#    #run transformix on points
#    points_file = point_transformix(pretransform_text_file, trfm_fl[-1], transformed_dst)
#    
#    #convert registered points into structure counts
#    converted_points = unpack_pnts(points_file, transformed_dst) 
#
def transformed_cells_to_ann(fld, ann, dst, fl_nm):
    """ consolidating to one function """
    
    dct = {}
    
    for fl in fld:
        converted_points = os.path.join(fl, "posttransformed_zyx_voxels.npy")
        if not os.path.exists(converted_points): 
            converted_points = os.path.join(fl, "transformed_points/posttransformed_zyx_voxels.npy")
        point_lst = transformed_pnts_to_allen_helper_func(np.load(converted_points), ann, order = "ZYX")
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
    df.to_pickle(os.path.join(dst, fl_nm))
    
    return os.path.join(dst, fl_nm)

#pma2aba_transformed = [os.path.join(atl_dst, xx) for xx in lr_brains]
trnsfrmd = [os.path.join(trnsfrm_dst, xx) for xx in lr_brains]

right = transformed_cells_to_ann(trnsfrmd, ann_right, dst, "right_side_no_prog_at_each_level_pma.p")
left = transformed_cells_to_ann(trnsfrmd, ann_left, dst, "left_side_no_prog_at_each_level_pma.p")

#import dict of cells by region
r_cells_regions = pckl.load(open(os.path.join(dst, "right_side_no_prog_at_each_level_pma.p"), "rb"), encoding = "latin1")
r_cells_regions = r_cells_regions.to_dict(orient = "dict")      

contra = {}; ipsi = {} #collect contra and ipsi frame
for k,v in r_cells_regions.items():
    if lr_dist[k] < 0:
        contra[k] = v
    else:
        ipsi[k] = v

#LEFT SIDE
l_cells_regions = pckl.load(open(os.path.join(dst, "left_side_no_prog_at_each_level_pma.p"), "rb"), encoding = "latin1")
l_cells_regions = l_cells_regions.to_dict(orient = "dict")      

for k,v in l_cells_regions.items():
    if lr_dist[k] > 0:
        contra[k] = v
    else:
        ipsi[k] = v
        
contra_df = pd.DataFrame(contra)
contra_df.to_csv(os.path.join(dst, "for_tp/nc_contra_counts_25_brains_pma.csv")) 

ipsi_df = pd.DataFrame(ipsi)
ipsi_df.to_csv(os.path.join(dst, "for_tp/nc_ipsi_counts_25_brains_pma.csv")) 

##########################################NO NEED TO RUN AGAIN IF ALREADY RUN ONCE################################################################
#%%

cells_regions_pth = os.path.join(dst, "for_tp/nc_contra_counts_25_brains_pma.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.020
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

def get_progeny(dic,parent_structure,progeny_list):
    """ 
    ---PURPOSE---
    Get a list of all progeny of a structure name.
    This is a recursive function which is why progeny_list is an
    argument and is not returned.
    ---INPUT---
    dic                  A dictionary representing the JSON file 191231_20180306_jg_bl6f_prv_14_488_049na_z7d5um_50msec_10povlp_16-41-17
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
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#then get only layers
layer56 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6a" or progen[-8:] == "layer 6b" or progen[-8:] == "Layer 6a" or progen[-8:] == "Layer 6b" or progen[-7:] == "layer 5" or progen[-7:] == "Layer 5":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer56.append(np.array(counts).sum(axis = 0))
layer56 = np.array(layer56)                

#then get only layers
layer23 = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-9:] == "layer 2/3" or progen[-9:] == "Layer 2/3":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer23.append(np.array(counts).sum(axis = 0))
layer23 = np.array(layer23)        

#calculate fraction of counts that come from layer 4,5,6, or 2/3
fracl56 = np.sum(layer56, axis = 0)/np.sum(counts_per_struct, axis = 0)
mean_fracl56_per_struct = np.nanmean(fracl56, axis = 0)

fracl23 = np.sum(layer23, axis = 0)/np.sum(counts_per_struct, axis = 0)
mean_fracl23_per_struct = np.nanmean(fracl23, axis = 0)

#%%

cells_regions_pth = os.path.join(dst, "for_tp/nc_ipsi_counts_25_brains_pma.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#then get only layers
layer56_ipsi = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6a" or progen[-8:] == "layer 6b" or progen[-8:] == "Layer 6a" or progen[-8:] == "Layer 6b" or progen[-7:] == "layer 5" or progen[-7:] == "Layer 5":
            counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    layer56_ipsi.append(np.array(counts).sum(axis = 0))
layer56_ipsi = np.array(layer56_ipsi)                

#%%
#get contra/ipsi ratios
layer56_contra = layer56

_ccontra = np.asarray([[np.sum(xx[:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in layer56_contra.T])
_cipsi = np.asarray([[np.sum(xx[:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in layer56_ipsi.T])
ratio = _ccontra/_cipsi
mean_ratio = np.mean(ratio, axis = 0)
std_ratio = np.std(ratio, axis = 0)

#separate by injection, vermis vs. hemisphere

func = lambda xx: 0 if xx < 3 else 1
#prv
primary_pool_vh = np.array([func(xx) for xx in primary_pool])
ratio_vermis = ratio[np.where(primary_pool_vh == 0)]
ratio_hem = ratio[np.where(primary_pool_vh == 1)]

mean_ratio_vermis = np.mean(ratio_vermis, axis = 0)
mean_ratio_hem = np.mean(ratio_hem, axis = 0)
std_ratio_vermis = np.std(ratio_vermis, axis = 0)
std_ratio_hem = np.std(ratio_hem, axis = 0)
#%%

#layer 5+6 p counts maps
layer56 = layer56
pcounts = np.array([xx/sum(xx) for xx in layer56.T])*100

#make % counts map like the h129 dataset (nc only for now)

## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (9,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})

#set colorbar features 
maxpcount = 30
whitetext = 3
label_coordsy, label_coordsx  = -0.37,0.5 #for placement of vertical labels
annotation_size = "x-small" #for the number annotations inside the heatmap
brain_lbl_size = "x-small"
yaxis = sois #for density by nc areas map

#sort inj fractions by primary lob
sort_pcounts = [pcounts[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_pcounts = np.array(list(itertools.chain.from_iterable(sort_pcounts)))
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")

ax = axes[1]
show = np.fliplr(sort_pcounts).T

vmin = 0
vmax = maxpcount
cmap = plt.cm.viridis
cmap.set_over("orange")
#colormap
bounds = np.linspace(vmin,vmax,((vmax-vmin)/5)+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("% of total neocortical counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="x-small")
ax.set_ylabel("Neocortical areas", fontsize="small")
ax.yaxis.set_label_coords(label_coordsy, label_coordsx)

ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=brain_lbl_size, ha="right")

plt.savefig(os.path.join(fig_dst, "pcounts_nc.pdf"), bbox_inches = "tight")

#%%

#make density map like the h129 dataset (nc only for now)
#get layer5/6 volumes
layer56_vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        if progen[-8:] == "layer 6a" or progen[-8:] == "layer 6b" or progen[-8:] == "Layer 6a" or progen[-8:] == "Layer 6b" or progen[-7:] == "layer 5" or progen[-7:] == "Layer 5":
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    layer56_vol.append(np.array(counts).sum(axis = 0))
layer56_vol = np.array(layer56_vol)        

density_l56 = np.array([xx/(layer56_vol[i]*(scale_factor**3)) for i, xx in enumerate(layer56)]).T

## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (10,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})

#set colorbar features 
maxdensity = 200
whitetext = 7
label_coordsy, label_coordsx  = -0.30,0.5 #for placement of vertical labels
annotation_size = "x-small" #for the number annotations inside the heatmap
brain_lbl_size = "small"
yaxis = sois #for density by nc areas map

#sort inj fractions by primary lob
sort_density = [density_l56[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over("darkred")
#colormap
norm = colors.PowerNorm(gamma=0.5)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")

ax = axes[1]
show = np.fliplr(sort_density).T

vmin = 0
vmax = maxdensity
cmap = plt.cm.viridis
cmap.set_over("orange")
#colormap
bounds = np.linspace(vmin,vmax,((vmax-vmin)/50)+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Density (Cells/$mm^3$)", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)
# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="x-small")
ax.set_ylabel("Neocortical areas", fontsize="small")
ax.yaxis.set_label_coords(label_coordsy, label_coordsx)

ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=brain_lbl_size, ha="right")

plt.savefig(os.path.join(fig_dst, "density_nc.pdf"), bbox_inches = "tight")

#%%
#glm?
#VARIABLES FOR GLM           
#POOL NC REGIONS!!

#make pcounts array
#total_counts_per_brain = np.sum(cell_counts_per_brain, axis=1)
#pcounts = np.asarray([(xx/total_counts_per_brain[i])*100 for i, xx in enumerate(cell_counts_per_brain)])

pcounts_pool = np.array([np.array([xx[0]+xx[1]+xx[2]+xx[4], xx[6], xx[5]+xx[7], 
                                          xx[8]+xx[9], xx[10], xx[11], xx[12], 
                                          xx[13]+xx[14], xx[15]+xx[16]]) for xx in pcounts])

#for display
regions = np.asarray(["Infralimbic, Prelimbic,\n Ant. Cingulate, Orbital",
       "Agranular insula", "Gustatory, Visceral", #removed frontal pole as it does not have any layer 5/6 neurons
       "Somatomotor, Somatosensory", "Retrosplenial", "Post. Parietal", "Visual",
       "Temporal, Auditory", "Peririhinal, Ectorhinal"])


X = np.array([brain/brain.sum() for brain in frac_of_inj_pool])
Y = pcounts_pool    
#%%

##  glm
c_mat = []
mat = []
pmat = []
mat_shuf = []
p_shuf = []
fit = []
fit_shuf = []

for itera in range(1000):
    if itera%100 == 0: print(itera)
    if itera == 0:
        shuffle = False
        inj = X.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        fit_shuf.append([])
        inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:]
    for count, region in zip(Y.T, regions):
        inj_ = inj[~np.isnan(count)]
        count = count[~np.isnan(count)]

        # intercept:
        inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
        
        glm = sm.GLM(count, inj_, family=sm.families.Poisson())
        res = glm.fit()
        
        coef = res.params[:-1]
        se = res.bse[:-1]
        pvals = res.pvalues[:-1] 

        val = coef/se

        if not shuffle:
            c_mat.append(coef)
            mat.append(val)
            pmat.append(pvals)
            fit.append(res.fittedvalues)
            
        elif shuffle:
            mat_shuf[-1].append(val)
            p_shuf[-1].append(pvals)
            fit_shuf[-1].append(res.fittedvalues)
        # inspect residuals
        """
        if not shuffle:
            plt.clf()
            plt.scatter(res.fittedvalues, res.resid)
            plt.hlines(0, res.fittedvalues.min(), res.fittedvalues.max())
            plt.title(region)
            plt.xlabel("Fit-vals")
            plt.ylabel("Residuals")
            plt.savefig(os.path.join(dst, "resid_inspection-{}.png").format(region))
        """
        
c_mat = np.array(c_mat)
mat = np.array(mat) # region x inj
mat_shuf = np.array(mat_shuf) # region x inj
pmat = np.array(pmat) # region x inj
p_shuf = np.array(p_shuf)
fit = np.array(fit)
fit_shuf = np.array(fit_shuf)

#%%
## display
fig = plt.figure(figsize=(6,5))
ax = fig.add_axes([.4,.1,.5,.8])

#set white text limit here
whitetext = 8
annotation_size = "medium"#annotation/number sizes

# map 1: weights
show = np.flipud(mat) # NOTE abs

vmin = 0
vmax = 10
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,((vmax-vmin)/2)+1)
#bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, 
                  boundaries=bounds, format="%d", shrink=0.3, aspect=10)
cb.set_label("Weight / SE", fontsize="small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

#annotations
for ri,row in enumerate(show): 
    for ci,col in enumerate(row):
        if col > whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize=annotation_size)
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize=annotation_size)   
            
# signif
sig = np.flipud(pmat) < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in np.flipud(regions)], fontsize="medium")

plt.savefig(os.path.join(fig_dst, "nc_glm.pdf"), bbox_inches = "tight")

#%%

#import matplotlib.ticker as ticker
#
###PRV TO HSV COMPARISON
### rank correlation plot
#
##import data
#main_data_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_contra_ipsi_counts_densities.p"
#data = pckl.load(open(main_data_pth, "rb"), encoding = "latin1")
#
#cell_counts_per_brain_left = data["cell_counts_per_brain_left"]
#cell_counts_per_brain_right = data["cell_counts_per_brain_right"]
#density_per_brain_left = data["density_per_brain_left"]
#density_per_brain_right = data["density_per_brain_right"]
#volume_per_brain_left = data["volume_per_brain_left"]
#lr_dist = data["lr_dist"]
#nc_areas = data["nc_areas"] #gives order of nc areas also
#scale_factor = 0.025
#
##preprocessing into contra/ipsi counts per brain, per structure
#scale_factor = 0.025
#nc_left_counts = cell_counts_per_brain_left
#nc_right_counts = cell_counts_per_brain_right
#nc_density_left = density_per_brain_left
#nc_density_right = density_per_brain_right
#
#lrv = list(lr_dist.values())
#lr_brains = list(lr_dist.keys())
#
##dct is just for my sanity, so im not mixing up brains
#_ccontra = []; _cipsi = []; _dcontra = []; _dipsi = []
#for i in range(len(lr_brains)):
#    if lrv[i] > 0: #right
#        #counts
#        _ccontra.append(nc_left_counts[i])
#        _cipsi.append(nc_right_counts[i])
#        #density
#        _dcontra.append(nc_density_left[i])
#        _dipsi.append(nc_density_right[i])
#    elif lrv[i] < 0: #left
#        #counts
#        _ccontra.append(nc_right_counts[i])
#        _cipsi.append(nc_left_counts[i])
#        #density
#        _dcontra.append(nc_density_right[i])
#        _dipsi.append(nc_density_left[i])
#
#_ccontra = np.asarray(_ccontra); _dcontra = np.asarray(_dcontra)
#_cipsi = np.asarray(_cipsi); _dipsi = np.asarray(_dipsi)
#
###optionally mask frontal areas
#nc_areas_mask = [False, False, False, False, False, True, True, True, True, True, True, True, True, True, True,True, True]
##h129 - ASCENDING ORDER
#hsv_counts = cell_counts_per_brain_left+cell_counts_per_brain_right
#hsv_pcounts = np.array([xx/sum(xx) for xx in hsv_counts])
#hsv_density = (hsv_counts)/((volume_per_brain_left*2) * (scale_factor ** 3))
##switched to descending order
#hsv_nc_areas_rank = np.array(nc_areas)[nc_areas_mask][np.argsort(np.median(hsv_density, axis = 0)[nc_areas_mask])][::-1]
#
##prv - ASCENDING ORDER
#prv_counts = cell_counts_per_brain
#prv_pcounts = np.array([xx/sum(xx) for xx in prv_counts])
#prv_density = density_per_brain
##switched to descending order
#prv_nc_areas_rank = np.array(nc_areas)[nc_areas_mask][np.argsort(np.median(prv_density, axis = 0)[nc_areas_mask])][::-1]
#
#hsv_ranks = [i+1 for i,nuc in enumerate(hsv_nc_areas_rank)]
#prv_ranks = []
##make proper ranks
#for soi in hsv_nc_areas_rank:
#    for i,s in enumerate(prv_nc_areas_rank):
#        if soi == s:
#            prv_ranks.append(i+1)
#
#X = hsv_ranks
#Y = prv_ranks
#
#results = sm.OLS(Y,sm.add_constant(X)).fit()
#
#mean_slope = results.params[1]
#mean_r2 = results.rsquared
#mean_intercept = results.params[0]
#
#fig = plt.figure(figsize=(10,5))
#ax = fig.add_axes([.4,.1,.5,.8])
#
##plot as scatter   
#ax.scatter(y = Y, x = X, s = 70)
#
##plot fit line
#ax.plot(mean_slope*range(len(X))+mean_intercept, '--k')    
#ytick_spacing = 1; xtick_spacing = 1
#ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
#
##plot nc area labels
#for i, txt in enumerate(hsv_nc_areas_rank):
#    ax.annotate(txt, (X[i], Y[i]), fontsize = "x-small")
#        
#ax.set_xlabel("H129 Density rank order")
#ax.set_ylabel("PRV Density rank order")
#
##make text box
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.3)
#
#textstr = "\n".join((
#    "slope: {:0.2f}".format(mean_slope),
#    "$R^2$: {:0.2f}".format(mean_r2)))
#
#ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#
#
#plt.savefig(os.path.join(dst, "prv_hsv_density_rank_order_median_no_frontal.pdf"), bbox_inches = "tight")
#
#plt.close()
#
##%%
#
##filter by vermis vs hemisphere
#
##h129
##injection site analysis
#inj_pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
#inj_dct = pckl.load(open(inj_pth, "rb"), encoding = "latin1")
#hsv_brains = inj_dct["brainnames"]
#hsv_primary_pool = inj_dct["primary_pool"]
#hsv_ak_pool = inj_dct["cb_regions_pool"]
#
#ak_vh = np.array(["Vermis", "Hemisphere"])
#func = lambda xx: 0 if xx < 3 else 1
#hsv_primary_pool_vh = np.array([func(xx) for xx in hsv_primary_pool])
#hsv_density_vermis = hsv_density[np.where(hsv_primary_pool_vh == 0)]
#hsv_density_hem = hsv_density[np.where(hsv_primary_pool_vh == 1)]
#
##prv
#prv_primary_pool_vh = np.array([func(xx) for xx in primary_pool])
#prv_density_vermis = density_per_brain[np.where(prv_primary_pool_vh == 0)]
#prv_density_hem = density_per_brain[np.where(prv_primary_pool_vh == 1)]
#
##VERMIS
#hsv_nc_areas_rank_vermis = np.array(nc_areas)[nc_areas_mask][np.argsort(np.median(hsv_density_hem, axis = 0)[nc_areas_mask])][::-1]
#prv_nc_areas_rank_vermis = np.array(nc_areas)[nc_areas_mask][np.argsort(np.median(prv_density_vermis, axis = 0)[nc_areas_mask])][::-1]
#
#hsv_ranks_vermis = [i+1 for i,nuc in enumerate(hsv_nc_areas_rank_vermis)]
#prv_ranks_vermis = []
##make proper ranks
#for soi in hsv_nc_areas_rank_vermis:
#    for i,s in enumerate(prv_nc_areas_rank_vermis):
#        if soi == s:
#            prv_ranks_vermis.append(i+1)
#
#
#X = hsv_ranks_vermis
#Y = prv_ranks_vermis
#
#results = sm.OLS(Y,sm.add_constant(X)).fit()
#
#mean_slope = results.params[1]
#mean_r2 = results.rsquared
#mean_intercept = results.params[0]
#
#fig = plt.figure(figsize=(10,5))
#ax = fig.add_axes([.4,.1,.5,.8])
#
##plot as scatter   
#ax.scatter(y = Y, x = X, s = 70)
#
##plot fit line
#ax.plot(mean_slope*range(len(X))+mean_intercept, '--k')    
#ytick_spacing = 1; xtick_spacing = 1
#ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
#
##plot nc area labels
#for i, txt in enumerate(hsv_nc_areas_rank_vermis):
#    ax.annotate(txt, (X[i], Y[i]), fontsize = "x-small")
#        
#ax.set_xlabel("H129 Density rank order (vermis)")
#ax.set_ylabel("PRV Density rank order (vermis)")
#
##make text box
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#
#textstr = "\n".join((
#    "slope: {:0.2f}".format(mean_slope),
#    "$R^2$: {:0.2f}".format(mean_r2)))
#
#ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#
#plt.savefig(os.path.join(dst, "prv_hsv_density_rank_order_vermis_median_no_frontal.pdf"), 
#            bbox_inches = "tight")
#
#plt.close()
#
##%%
#
#import scipy, seaborn as sns, json
#
#    # Now write the function to get all progeny of an input nodename
#def get_progeny(dic,parent_structure,progeny_list):
#    """ 
#    ---PURPOSE---
#    Get a list of all progeny of a structure name.
#    This is a recursive function which is why progeny_list is an
#    argument and is not returned.
#    ---INPUT---
#    dic                  A dictionary representing the JSON file 
#                         which contains the ontology of interest
#    parent_structure     The structure
#    progeny_list         The list to which this function will 
#                         append the progeny structures. 
#    """
#    if 'msg' in list(dic.keys()): dic = dic['msg'][0]
#    
#    name = dic.get('name')
#    children = dic.get('children')
#    if name == parent_structure:
#        for child in children: # child is a dict
#            child_name = child.get('name')
#            progeny_list.append(child_name)
#            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
#        return
#    
#    for child in children:
#        child_name = child.get('name')
#        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
#    return 
#
##get progeny of all large structures
#ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"
#
#with open(ontology_file) as json_file:
#    ontology_dict = json.load(json_file)
#
##rank by cfos mean cells /mm3 vs mean cells /mm3. Or make cfos/mm3 delta between the two..
#cs = pd.read_pickle("/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2018_sfn_poster/cfos/stat_dataframe.p")
#
##adding up to hierarchy
##add parent regions to df
#df = pd.DataFrame(columns = cs.columns)
#df["Structure"] = nc_areas
#cs = cs.append(df, ignore_index = True)
#
#for soi in nc_areas:
#    progeny = []
#    get_progeny(ontology_dict, soi, progeny)
#    counts = [cs.loc[cs.Structure == xx, "Stimulation mean"].values[0] for xx in progeny if len(cs.loc[cs.Structure == xx, "Stimulation mean"].values)]
#    vol = [cs.loc[cs.Structure == xx, "structure_vol"].values[0] for xx in progeny if len(cs.loc[cs.Structure == xx, "structure_vol"].values)]
#    counts = sum(counts); vol = sum(vol)
#    cs.loc[cs.Structure == soi, "Stimulation mean"] = counts
#    cs.loc[cs.Structure == soi, "structure_vol"] = vol
#
#cfos = cs[['Structure','Stimulation mean', "structure_vol"]].copy()
#
#cfos["Stimulation density mean"] = cfos["Stimulation mean"]/(cfos["structure_vol"]/(0.020**3))
#
#cfos = cfos[cfos.Structure.isin(nc_areas)]
#cfos = cfos.sort_values('Stimulation density mean', ascending=False)
#cfos['cfos_order'] = range(1,len(cfos)+1)
#
#prv_nc_areas_rank = np.array(nc_areas)[np.argsort(np.median(prv_density, axis = 0))][::-1]
#cfos_nc_areas_rank = cfos.Structure.values
##%%
#
#cfos_ranks = [i+1 for i,nuc in enumerate(cfos_nc_areas_rank)]
#prv_ranks = []
##make proper ranks
#for soi in cfos_nc_areas_rank :
#    for i,s in enumerate(prv_nc_areas_rank):
#        if soi == s:
#            prv_ranks.append(i+1)
#
#
#X = prv_ranks
#Y = cfos_ranks
#
#df = pd.DataFrame()
#df["cfos_ranks"] = cfos_ranks
#df["prv_ranks"] = prv_ranks
#
#slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df[['prv_ranks', 'cfos_ranks']].values)
#
#g = sns.regplot(data=df, x='prv_ranks', y = 'cfos_ranks')
#g.text(19,17,'slope {:0.2f}\nr_value {:0.2f}\np_value {:0.2f}'.format(slope, r_value, p_value))
#
#ytick_spacing = 1; xtick_spacing = 1
#g.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
#g.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
#
##plot nc area labels
#for i, txt in enumerate(cfos_nc_areas_rank):
#    g.annotate(txt, (X[i], Y[i]), fontsize = "x-small")
#
#plt.savefig(os.path.join(dst, 'cfos_v_prv.jpg'), bbox_inches = "tight", dpi=300); plt.close()
#with open(os.path.join(dst, 'all_nc_structures_cfosbyactivationratio.txt'), 'a') as fl:
#    fl.write('slope {}, intercept {},\n r_value {}, p_value {},'.format(slope, intercept, r_value, p_value))
#    fl.close()