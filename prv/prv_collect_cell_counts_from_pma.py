#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 15:44:38 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl, statsmodels.api as sm, itertools
from skimage.external import tifffile
from scipy.ndimage.measurements import center_of_mass

from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister, change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files, point_transformix, unpack_pnts
from tools.utils.io import makedir
from tools.analysis.network_analysis import make_structure_objects

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#set paths, variables, etc.
dst = "/jukebox/wang/zahra/prv/"
pma_ann_pth = os.path.join(dst, "pma_annotation_sagittal_atlas_20um_iso_60um_edge_160um_vntric_erosion.tif")
pma_ann = tifffile.imread(pma_ann_pth)

#collect 
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20180205_jg_bl6f_prv_02", "20180205_jg_bl6f_prv_03", "20180205_jg_bl6f_prv_04", "20180215_jg_bl6f_prv_05", "20180215_jg_bl6f_prv_06",
       "20180215_jg_bl6f_prv_09", "20180305_jg_bl6f_prv_12", "20180305_jg_bl6f_prv_15", "20180312_jg_bl6f_prv_17", "20180326_jg_bl6f_prv_37",
       "20180313_jg_bl6f_prv_21", "20180313_jg_bl6f_prv_23", "20180313_jg_bl6f_prv_24", "20180305_jg_bl6f_prv_11", "20180313_jg_bl6f_prv_25",
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
ak_pool = np.array(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X",
                 "Simplex", "Crura", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]),xx[10],xx[11:13].sum(),np.sum(xx[13:16])] for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
#get n"s after pooling
primary_lob_n = np.array([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])


#make structures
#FIXME: for some reason the allen table does not work on this, is it ok to use PMA        
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=pma_ann_pth)

#set variables
lr_brains = list(lr_dist.keys())

atl_dst = os.path.join(dst, "pma_to_aba"); makedir(atl_dst)
trnsfrm_dst = os.path.join(dst, "prv_transformed_cells")
id_table = pd.read_excel(df_pth)

#%%
cells_src = os.path.join(dst, "prv_transformed_cells")
post_transformed = [os.path.join(cells_src, os.path.join(xx, "transformed_points/posttransformed_zyx_voxels.npy")) for xx in lr_brains]
transformfiles = ["/jukebox/wang/zahra/aba_to_pma/TransformParameters.0.txt",
                  "/jukebox/wang/zahra/aba_to_pma/TransformParameters.1.txt"]
##########################################NO NEED TO RUN AGAIN IF ALREADY RUN ONCE################################################################
#collect 
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

whl_brain = transformed_cells_to_ann(trnsfrmd, pma_ann, dst, "whl_brain_no_prog_at_each_level_pma_atl.p")
##########################################NO NEED TO RUN AGAIN IF ALREADY RUN ONCE################################################################
#%%
def get_cell_n_density_counts(brains, structure, structures, cells_regions, id_table, scale_factor = 0.025):
    """ consolidating to one function bc then no need to copy/paste """
    #get cell counts adn densities
    #get densities for all the structures
    df = pd.read_excel(id_table, index_col = None)
    df = df.drop(columns = ["Unnamed: 0"])
    df = df.sort_values(by = ["name"])
    
    #make new dict - for all brains
    cells_pooled_regions = {} #for raw counts
    volume_pooled_regions = {} #for density
    
    for brain in brains:    
        #make new dict - this is for EACH BRAIN
        c_pooled_regions = {}
        d_pooled_regions = {}
        
        for soi in structure:
            print(soi)
            try:
                soi = [s for s in structures if s.name==soi][0]
                counts = [] #store counts in this list
                volume = [] #store volume in this list
                for k, v in cells_regions[brain].items():
                    if k == soi.name:
                        counts.append(v)
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))
                progeny = [str(xx.name) for xx in soi.progeny]
                #now sum up progeny
                if len(progeny) > 0:
                    for progen in progeny:
                        for k, v in cells_regions[brain].items():
                            if k == progen and progen != "Primary somatosensory area, unassigned, layer 4,5,6":
                                counts.append(v)
                                #add to volume list from LUT
                                volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0])
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
            except Exception as e:
                print(e)
                for k, v in cells_regions[brain].items():
                    if k == soi:
                        counts.append(v)                    
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
                        
        #add to big dict
        cells_pooled_regions[brain] = c_pooled_regions
        volume_pooled_regions[brain] = d_pooled_regions
    #making the proper array per brain where regions are removed
    cell_counts_per_brain = []
    #initialise dummy var
    i = []
    for k,v in cells_pooled_regions.items():
        dct = cells_pooled_regions[k]
        for j,l in dct.items():
            i.append(l)  
        cell_counts_per_brain.append(i)
        #re-initialise for next
        i = []  
    cell_counts_per_brain = np.asarray(cell_counts_per_brain)
    
    volume_per_brain = []
    #initialise dummy var
    i = []
    for k,v in volume_pooled_regions.items():
        dct = volume_pooled_regions[k]
        for j,l in dct.items():
            i.append(l)  
        volume_per_brain.append(i)
        #re-initialise for next
        i = []  
    volume_per_brain = np.asarray(volume_per_brain)
    #calculate denisty
    density_per_brain = np.asarray([xx/(volume_per_brain[i]*(scale_factor**3)) for i, xx in enumerate(cell_counts_per_brain)])
    
    return cell_counts_per_brain, density_per_brain, volume_per_brain

#making dictionary of cells by region
cells_regions = pckl.load(open(os.path.join(dst, "whl_brain_no_prog_at_each_level_pma_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

areas = ['Paraventricular hypothalamic nucleus', 'Periventricular hypothalamic nucleus, posterior part',
       'Periventricular hypothalamic nucleus, anterior part', 'Periventricular hypothalamic nucleus, intermediate part',
       'Periventricular hypothalamic nucleus, preoptic part', 'Arcuate hypothalamic nucleus', 
       'Anterodorsal preoptic nucleus', 'Anteroventral preoptic nucleus',
       'Anteroventral periventricular nucleus',
       'Dorsomedial nucleus of the hypothalamus',
       'Median preoptic nucleus', 'Medial preoptic area', 'Ventrolateral preoptic nucleus',
       'Anterior hypothalamic nucleus', 'Mammillary body',
       'Medial mammillary nucleus',
       'Supramammillary nucleus', 'Tuberomammillary nucleus',
       'Medial preoptic nucleus', 'Dorsal premammillary nucleus',
       'Ventral premammillary nucleus', 'Ventromedial hypothalamic nucleus', 'Posterior hypothalamic nucleus',
       'Lateral hypothalamic area', 'Lateral preoptic area', 'Preparasubthalamic nucleus', 'Parasubthalamic nucleus',
       'Retrochiasmatic area','Subthalamic nucleus', 'Tuberal nucleus', 'Zona incerta', 'Median eminence']

#RIGHT SIDE
cell_counts_per_brain, density_per_brain, volume_per_brain = get_cell_n_density_counts(brains, areas, structures, 
                                                                                       cells_regions, df_pth, scale_factor = 0.020)

#%%
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (10,10), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1,5]})

#set colorbar features 
maxdensity = 1500
label_coordsy, label_coordsx  = -0.50,0.5 #for placement of vertical labels
annotation_size = "x-small" #for the number annotations inside the heatmap
brain_lbl_size = "small"
yaxis = areas #for density by nc areas map
    
#sort inj fractions by primary lob
sort_density = [density_per_brain[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
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
cmap.set_over('darkred')
#colormap
bounds = np.linspace(vmin,vmax,6)
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
show = np.fliplr(sort_density).T

vmin = 0
vmax = maxdensity
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.2, aspect=10)
cb.set_label("Density (Cells/$mm^3$)", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="x-small")
ax.set_ylabel("Hypothalamic areas", fontsize="small")
ax.yaxis.set_label_coords(label_coordsy, label_coordsx)

ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=brain_lbl_size, ha="right")

plt.savefig(os.path.join(dst, "density_hypothal.pdf"), bbox_inches = "tight")

#%%
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (10,10), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1,5]})

#set colorbar features 
maxpcounts = 10
label_coordsy, label_coordsx  = -0.50,0.5 #for placement of vertical labels
annotation_size = "x-small" #for the number annotations inside the heatmap
brain_lbl_size = "small"
yaxis = areas #for density by nc areas map

#make pcounts array
total_counts_per_brain = np.sum(cell_counts_per_brain, axis=1)
pcounts = np.asarray([(xx/total_counts_per_brain[i])*100 for i, xx in enumerate(cell_counts_per_brain)])
#sort inj fractions by primary lob
sort_pcounts = [pcounts[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_pcounts = np.array(list(itertools.chain.from_iterable(sort_pcounts)))
    
#sort inj fractions by primary lob
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
bounds = np.linspace(vmin,vmax,6)
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
vmax = maxpcounts
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.2, aspect=10)
cb.set_label("Density (Cells/$mm^3$)", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="x-small")
ax.set_ylabel("Hypothalamic areas", fontsize="small")
ax.yaxis.set_label_coords(label_coordsy, label_coordsx)

ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=brain_lbl_size, ha="right")

plt.savefig(os.path.join(dst, "pcounts_hypothal.pdf"), bbox_inches = "tight")
