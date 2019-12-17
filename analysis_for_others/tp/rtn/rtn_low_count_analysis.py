#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 11:14:42 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl, itertools
from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.registration.register import change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files
from tools.registration.transform_list_of_points import point_transformix, unpack_pnts
from tools.utils.io import makedir, load_kwargs
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass
import matplotlib.colors, statsmodels.api as sm
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop
from scipy.stats import median_absolute_deviation as mad
from scipy.ndimage.filters import gaussian_filter as gfilt
from scipy.ndimage import label

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

def find_site(im, thresh=3, filter_kernel=(3,3,3), num_sites_to_keep=1):

    if type(im) == str: im = tifffile.imread(im)

    filtered = gfilt(im, filter_kernel)
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)

    if nlab == 0:
        raise Exception('Site not detected, try a lower threshold?')
    elif nlab == 1:
        return labelled.astype(bool)
    elif num_sites_to_keep == 1:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        return labelled == np.argmax(sizes)+1
    else:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
        return np.in1d(labelled, vals).reshape(labelled.shape)
    
dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/rtn"
fig_dst = "/home/wanglab/Desktop"
ann_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/atlases/sagittal_allen_ann_25um_iso_60um_edge_160um_ventricular_erosion.tif"

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
brains = ["20160622_db_bl6_crii_52hr_01", "20160622_db_bl6_unk_01", "20160801_db_cri_02_1200rlow_52hr", 
          "20160801_db_l7_cri_01_mid_64hr", "20160822_tp_bl6_crii_250r_01", "20160822_tp_bl6_crii_1250r_05", 
          "20160822_tp_bl6_crii_1500r_06", "20160823_tp_bl6_cri_250r_01",
          "20160823_tp_bl6_cri_500r_02", "20160912_tp_bl6_lob7_750r_04", "20160916_tp_lob6_ml_01", "20160916_tp_lob6_ml_04", 
          "20160916_tp_lob7_250r_05", "20160920_tp_bl6_lob7_250r_02", "20160920_tp_bl6_lob7_500r_03", "20160920_tp_bl6_lob7_ml_01", 
          "20161201_db_bl6_lob6b_500r_53d5hr", "20161203_tp_bl6_crii_250r_08", 
          "20161203_tp_bl6_lob7_500r_05", "20161203_tp_bl6_lob7_1000r_06", "20161203_tp_bl6_lob7_1500r_07",
          "20161203_tp_bl6_lob7_ml_04", "20161205_tp_bl6_sim_250r_02", "20161205_tp_bl6_sim_750r_03",
          "20161205_tp_bl6_sim_1250r_04", "20161205_tp_bl6_sim_1750r_05","20161207_db_bl6_lob6a_500r_53hr", 
          "20161207_db_bl6_lob6a_850r_53hr", "20161208_db_bl6_cri_50r_pv_53hr","20170115_tp_bl6_lob6b_500r_05",
          "20170116_tp_bl6_lob6b_lpv_07","20170130_tp_bl6_sim_rlat_05",
          "20170207_db_bl6_crii_1300r_02", "20170212_tp_bl6_lob8_ml_05",
          "20170308_tp_bl6f_cri_2x_03","20170308_tp_bl6f_lob7_2x_02", 
          "20170410_tp_bl6_lob6a_ml_repro_01", "20170410_tp_bl6_lob6a_ml_repro_02",
          "20170411_db_bl6_crii_lat_53hr", "20170411_db_bl6_crii_mid_53hr",
          "20170411_db_bl6_crii_rpv_53hr", "20170419_db_lob6b_rpv_53hr", "20170419_db_bl6_cri_rpv_53hr", "20170207_db_bl6_crii_rpv_01",
        "20170419_db_bl6_cri_mid_53hr", "20161207_db_bl6_lob6a_50rml_53d5hr", "20161205_tp_bl6_lob45_1000r_01"]    

brains_w_no_inj_vol = ["20170308_tp_bl6f_cri_2x_03", "20160920_tp_bl6_lob7_ml_01", "20170207_db_bl6_crii_1300r_02",
                       "20160801_db_cri_02_1200rlow_52hr", "20170115_tp_bl6_lob6b_500r_05"]

src = "/jukebox/wang/pisano/tracing_output/antero_4x"

imgs = [os.path.join(src, xx) for xx in brains]

#pool brain names and L/R designation into dict
lr_dist = {}
thal_inj_vol = {}
brains_w_no_inj = []

#save inj segments in folder
sv_dst = os.path.join(dst, "injection"); makedir(sv_dst)
#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    kwargs = load_kwargs(img)
    
    print(brain)
    try:
        if brain in brains_w_no_inj_vol:
            inj_vol_pth = os.path.join(dst, brain+"_reg/result.1.tif")
            inj_vol = tifffile.imread(inj_vol_pth)
        else:
            inj_vol_pth = [xx for xx in kwargs["volumes"] if xx.ch_type == "injch"][0]
            inj_vol = tifffile.imread(inj_vol_pth.ch_to_reg_to_atlas)
        inj_vol[inj_vol < 0] = 10000
        assert np.sum(inj_vol < 0) == 0
        z,y,x = inj_vol.shape
        
        arr = find_site(inj_vol[:, 423:, :])
        #save segment
        tifffile.imsave(os.path.join(sv_dst, brain+".tif"), arr.astype("uint16"))
        
        z_c,y_c,x_c = center_of_mass(arr)
        #take distance from center to arbitrary "midline" (aka half of z axis)
        dist = z_c-(z/2)
        #save to dict 
        lr_dist[brain] = dist
        thal_inj_vol[brain] = np.sum(inj_vol)
        
        if dist < 0:
            print("brain {} has a left-sided injection\n".format(brain))
        elif dist > 0:
            print("brain {} has a right-sided injection\n".format(brain))
        else:
            print("brain has an injection close to midline so not considering it rn\n")
    except:
        print("brain %s has no injection volume, segment from elsewhere\n" % brain)
        brains_w_no_inj.append(img)

#FIXME: temporarily dropping these brains with no inj, might use them again later
#%%
#make structures
#FIXME: for some reason the allen table does not work on this, is it ok to use PMA..?
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"

structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
lr_brains = list(lr_dist.keys())
atl_dst = os.path.join(dst, "pma_to_aba"); makedir(atl_dst)
id_table = pd.read_excel(df_pth)

#%%
#------------------------------------------------------------------------------------------------------------------------------
#NOTE THAT ONLY HAVE TO DO THIS ONCE!!!! DO NOT NEED TO DO AGAIN UNLESS DOUBLE CHECKIHG
#transform points to allen atlas space

#get brains that we actually need to get cell counts from
#src = "/jukebox/scratch/zmd/h129_qc/rtn_thal_transformed_points"
#post_transformed = [os.path.join(src, os.path.join(xx, "transformed_points/posttransformed_zyx_voxels.npy")) for xx in lr_brains]
#transformfiles = ["/jukebox/wang/zahra/aba_to_pma/TransformParameters.0.txt",
#                  "/jukebox/wang/zahra/aba_to_pma/TransformParameters.1.txt"]
#
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
    
#%%
#------------------------------------------------------------------------------------------------------------------------------    
#NOTE THAT ONLY HAVE TO DO THIS ONCE!!!! DO NOT NEED TO DO AGAIN UNLESS DOUBLE CHECKIHG
#def transformed_cells_to_ann(fld, ann, dst, fl_nm):
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
#right = transformed_cells_to_ann(pma2aba_transformed, ann_right, dst, "thal_rtn_right_side_no_prog_at_each_level_allen_atl.p")
##collect counts from left side
#left = transformed_cells_to_ann(pma2aba_transformed, ann_left, dst, "thal_rtn_left_side_no_prog_at_each_level_allen_atl.p")

#%%
def get_cell_n_density_counts(brains, structure, structures, cells_regions, scale_factor = 0.025):
    """ consolidating to one function bc then no need to copy/paste """
    #get cell counts adn densities
    #get densities for all the structures
    df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts_16bit.xlsx", index_col = None)
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
            try:
                soi = [s for s in structures if s.name==soi][0]
                counts = [] #store counts in this list
                volume = [] #store volume in this list
                for k, v in cells_regions[brain].items():
                    if k == soi.name:
                        counts.append(v)
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0]/2)#divide by 2 bc its half brain
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
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0]/2)#divide by 2 bc its half brain
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
cells_regions = pckl.load(open(os.path.join(dst, "thal_rtn_right_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

nuclei = ["Thalamus", "Ventral posteromedial nucleus of the thalamus", "Ventral posterolateral nucleus of the thalamus",
          "Ventral anterior-lateral complex of the thalamus", "Ventral medial nucleus of the thalamus", "Anteroventral nucleus of thalamus", 
          "Reticular nucleus of the thalamus", "Ventral part of the lateral geniculate complex", "Mediodorsal nucleus of thalamus",
          "Submedial nucleus of the thalamus", "Nucleus of reuniens", "Paraventricular nucleus of the thalamus", 
          "Central lateral nucleus of the thalamus", "Parafascicular nucleus", "Posterior complex of the thalamus",
          "Lateral dorsal nucleus of thalamus", "Lateral posterior nucleus of the thalamus", "Lateral habenula"]

#RIGHT SIDE
cell_counts_per_brain_right, density_per_brain_right, volume_per_brain_right = get_cell_n_density_counts(lr_brains, 
                                                                               nuclei, structures, cells_regions)
#LEFT SIDE
cells_regions = pckl.load(open(os.path.join(dst, "thal_rtn_left_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      
cell_counts_per_brain_left, density_per_brain_left, volume_per_brain_left = get_cell_n_density_counts(lr_brains, 
                                                                               nuclei, structures, cells_regions)


#regression total count on reticular thalamus count (y) for these brains
#will sum left and right for this
#plt.figure()
##X = np.sum(cell_counts_per_brain_left, axis=1)+np.sum(cell_counts_per_brain_right, axis=1)
#x = np.sum(cell_counts_per_brain_left, axis=1)+np.sum(cell_counts_per_brain_right, axis=1)
##mask high count brain
#X = x[x < 1000]
##plot vpm and vpl alongside also
#lbls = ["VPL", "VPM", "RTN"]
#cols = ["r", "g", "b"]
#idxs = [2, 3]
#for i,idx in enumerate(idxs):
#    Y = cell_counts_per_brain_right[:,idx]+cell_counts_per_brain_left[:,idx]
#    #mask 
#    Y = Y[x < 1000]
#    plt.scatter(x = X, y = Y, facecolors="none", edgecolors=cols[i], label = lbls[i])
#
#plt.xlabel("Total thalamic counts")    
#plt.ylabel("Thalamic nuclei counts")
#plt.xlim([0, 200])
#plt.ylim([0, 30])
#plt.legend()
#plt.savefig(os.path.join(fig_dst, "thal_nuclei_reg.pdf"))


#get inj fractions, using GLM code
#making dictionary of injection sites
injections = {}

for i,pth in enumerate(lr_brains):
    print(i, pth)
    pth = os.path.join(sv_dst, pth+".tif")
    injection = tifffile.imread(pth)
    print(injection.shape)
    injections[os.path.basename(pth)] = injection #files have 2 .tif in the end
    
inj_raw = np.array([inj.astype(bool) for nm, inj in injections.items()])
    
#use original annotation file
ann_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")[:, 423:, :] #based on seg script
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)

#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lobule III": 984,
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

#get ns
primary_lob_n = np.asarray([np.where(primary == i)[0].shape[0] for i in np.unique(primary)])

#%%

#make injection heatmap figure
from tools.imageprocessing.orientation import fix_orientation

atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"

sites = np.array([fix_orientation(xx, ("2","0","1")) for i,xx in enumerate(inj_raw) 
                    if not lr_brains[i] == "20170419_db_bl6_cri_mid_53hr"]) #the y-axis cutoff for visualization
atl = fix_orientation(tifffile.imread(atl_pth)[:, 450:, :], ("2","0","1"))

my_cmap = eval("plt.cm.{}(np.arange(plt.cm.RdBu.N))".format("viridis"))
my_cmap[:1,:4] = 0.0  
my_cmap = mpl.colors.ListedColormap(my_cmap)
my_cmap.set_under("w")
plt.figure()
plt.imshow(np.max(atl, axis=0), cmap="gray")
plt.imshow(np.max(np.sum(sites, axis=0), axis = 0), alpha=0.90, cmap=my_cmap); plt.colorbar(); plt.axis("off")

plt.savefig(os.path.join(fig_dst, "inj_heatmap.pdf"), dpi = 300, transparent = True);
plt.close()
#%%
#regression total count on reticular thalamus count (y) for these brains
#will sum left and right for this
#ONLY VERMIS
#X = (np.sum(cell_counts_per_brain_left, axis=1)+np.sum(cell_counts_per_brain_right, axis=1))[primary < 10]
#plt.figure()
#x = (np.sum(cell_counts_per_brain_left, axis=1)+np.sum(cell_counts_per_brain_right, axis=1))[primary < 10]
##mask high count brain
#X = x[x < 1000]
##plot vpm and vpl alongside also
#lbls = ["VPL", "VPM", "RTN"]
#cols = ["r", "g", "b"]
#idxs = [2, 3, 13]
#for i,idx in enumerate(idxs):
#    Y = (cell_counts_per_brain_right[:,idx]+cell_counts_per_brain_left[:,idx])[primary < 10]
#    #mask 
#    Y = Y[x < 1000]
#    plt.scatter(x = X, y = Y, facecolors="none", edgecolors=cols[i], label = lbls[i])
#
#plt.xlabel("Total thalamic counts")    
#plt.ylabel("Thalamic nuclei counts")
#plt.title("Only vermis injections")
##plt.xlim([0, 200])
##plt.ylim([0, 30])
#plt.legend()
#plt.savefig(os.path.join(fig_dst, "thal_nuclei_reg_only_vermis.pdf"))
#
#
##%%
##get ratios of vpm/rtn over total thalamic counts
#fig = plt.figure(figsize=(8,5))
#ax = fig.add_axes([.4,.1,.5,.8])
#
#X = (np.sum(cell_counts_per_brain_left, axis=1)+np.sum(cell_counts_per_brain_right, axis=1))[cell_counts_per_brain_left[:,13] > 0]
#Y = (cell_counts_per_brain_right[:,2:4].sum(axis=1)/cell_counts_per_brain_left[:,13])[cell_counts_per_brain_left[:,13] > 0]
#results = sm.OLS(Y,sm.add_constant(X)).fit()
#
#plt.scatter(x = X, y = Y, facecolors="none", edgecolors="k")
#plt.xlabel("Total thalamic counts")    
#plt.ylabel("VPM+VPL/RTN counts")
##plt.xlim([0, 200])
##plt.ylim([0, 30])
#
#textstr = "\n".join((
#    "slope: {:0.5f}".format(results.params[1]),
#    "$R^2$: {:0.5f}".format(results.rsquared)))
#
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
## place a text box in upper left in axes coords
#ax.text(0.7, 0.95, textstr, transform=ax.transAxes, fontsize=12,
#        verticalalignment='top', bbox=props)
##plt.xlim([0, 400])
##plt.ylim([0, 20])
#
#plt.savefig(os.path.join(fig_dst, "thal_rtn_ratio.pdf"))
#%%
#preprocessing into contra/ipsi counts per brain, per structure
scale_factor = 0.025
nc_left_counts = cell_counts_per_brain_left
nc_right_counts = cell_counts_per_brain_right
nc_density_left = density_per_brain_left
nc_density_right = density_per_brain_right

lrv = list(lr_dist.values())
lr_brains = list(lr_dist.keys())

#dct is just for my sanity, so im not mixing up brains
_ccontra = []; _cipsi = []; _dcontra = []; _dipsi = []
for i in range(len(lr_brains)):
    if lrv[i] > 0: #right
        #counts
        _ccontra.append(nc_left_counts[i])
        _cipsi.append(nc_right_counts[i])
        #density
        _dcontra.append(nc_density_left[i])
        _dipsi.append(nc_density_right[i])
    elif lrv[i] < 0: #left
        #counts
        _ccontra.append(nc_right_counts[i])
        _cipsi.append(nc_left_counts[i])
        #density
        _dcontra.append(nc_density_right[i])
        _dipsi.append(nc_density_left[i])


_ccontra = np.asarray(_ccontra).T; _dcontra = np.asarray(_dcontra).T
_cipsi = np.asarray(_cipsi).T; _dipsi = np.asarray(_dipsi).T
_ratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
#make into one
_dist = np.asarray(list(lr_dist.values()))


_inj = expr_all_as_frac_of_inj
 
_primary = np.asarray([primary[i] for i in range(len(primary)) if brains[i] in lr_brains])

#sort by distance
sort_dist = np.sort(_dist)
sort_ccontra = _ccontra.T[np.argsort(_dist, axis = 0)]
sort_cipsi = _cipsi.T[np.argsort(_dist, axis = 0)]
sort_dcontra = _dcontra.T[np.argsort(_dist, axis = 0)]
sort_dipsi = _dipsi.T[np.argsort(_dist, axis = 0)]
sort_ratio = _ratio.T[np.argsort(_dist, axis = 0)]
sort_vox_per_region = volume_per_brain_left[np.argsort(_dist, axis = 0)]
sort_inj = _inj[np.argsort(_dist)]   
sort_brains = np.array(lr_brains)[np.argsort(_dist)]

print(sort_dist.shape)
print(sort_ratio.shape)

#group thalamus regions in smaller, meta regions
grps = np.array(["Sensory-motor thalamus" , "Polymodal thalamus"])
sort_ccontra_pool = np.asarray([[np.sum(xx[:8]), np.sum(xx[8:-1])] for xx in sort_ccontra])
sort_dcontra_pool = np.asarray([[np.sum(xx[:8]), np.sum(xx[8:-1])] for xx in sort_ccontra])/(np.asarray([[np.sum(xx[:8]), 
                                 np.sum(xx[8:])] for xx in sort_vox_per_region])*(scale_factor**3))
sort_cipsi_pool = np.asarray([[np.sum(xx[:8]), np.sum(xx[8:-1])] for xx in sort_cipsi])
sort_dipsi_pool = np.asarray([[np.sum(xx[:8]), np.sum(xx[8:-1])] for xx in sort_cipsi])/(np.asarray([[np.sum(xx[:8]), 
                                 np.sum(xx[8:])] for xx in sort_vox_per_region])*(scale_factor**3))
sort_cratio_pool = np.asarray([sort_ccontra_pool[i]/sort_cipsi_pool[i] for i in range(len(sort_ccontra_pool))])
sort_dratio_pool = np.asarray([sort_dcontra_pool[i]/sort_dipsi_pool[i] for i in range(len(sort_dcontra_pool))])

#pooled injections
ak_pool = np.asarray(["Lob. III, IV-V", "Lob. VIa, VIb, VII-X", 
                 "Simplex", "Crus I", "Crus II", "PM, CP"])

#pooling injection regions
expr_all_as_frac_of_lob_pool = np.asarray([[xx[0]+xx[1], xx[2]+xx[3]+xx[4]+xx[5]+xx[6]+xx[7], 
                                            xx[8], xx[9], xx[10], xx[11]+xx[12]] for xx in expr_all_as_frac_of_lob])

expr_all_as_frac_of_inj_pool = np.asarray([[xx[0]+xx[1], xx[2]+xx[3]+xx[4]+xx[5]+xx[6]+xx[7], 
                                            xx[8], xx[9], xx[10], xx[11]+xx[12]] for xx in expr_all_as_frac_of_inj])
    
primary_pool = np.asarray([np.argmax(e) for e in expr_all_as_frac_of_inj_pool])
primary_lob_n = np.asarray([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])

sort_inj_pool = np.array([expr_all_as_frac_of_lob_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)])
sort_inj_pool = np.array(list(itertools.chain.from_iterable(sort_inj_pool)))

#%%

#make % counts map like the h129 dataset (nc only for now)

#drop brains which have < 10 cells in thalamus
lowest_count = 10
mask = np.array([xx[0] for xx in cell_counts_per_brain_left[:, :1] > lowest_count])
sort_inj_pool_filtered = sort_inj_pool[mask]
primary_pool_filtered = primary_pool[mask]

## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (15,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})

#set colorbar features 
maxpcount = 30
whitetext = 3
label_coordsy, label_coordsx  = -0.37,0.5 #for placement of vertical labels
annotation_size = "x-small" #for the number annotations inside the heatmap
brain_lbl_size = "x-small"
yaxis = nuclei[1:] #for density by nc areas map

#make pcounts array
cell_counts_per_brain = (cell_counts_per_brain_left+cell_counts_per_brain_right)[mask]
total_counts_per_brain = cell_counts_per_brain_left[:, :1][mask]
pcounts = np.asarray([(xx/total_counts_per_brain[i])*100 for i, xx in enumerate(cell_counts_per_brain[:, 1:])])
#sort inj fractions by primary lob
sort_pcounts = [pcounts[np.where(primary_pool_filtered == idx)[0]] for idx in np.unique(primary_pool_filtered)]
sort_pcounts = np.array(list(itertools.chain.from_iterable(sort_pcounts)))

sort_brains = [np.array(lr_brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = np.array(list(itertools.chain.from_iterable(sort_brains)))[mask]

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj_pool_filtered).T

vmin = 0
vmax = 0.3
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
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,((vmax-vmin)/5)+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("% of total thalamic counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="x-small")
ax.set_ylabel("Thalamic nuclei", fontsize="small")
ax.yaxis.set_label_coords(label_coordsy, label_coordsx)

ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=brain_lbl_size, ha="right")

plt.savefig(os.path.join(fig_dst, "pcounts_thal.pdf"), bbox_inches = "tight")

#%%
#filter results by contra/ipsi ratio > 1 in thalamus

threshold = 1

#drop brains which have < 10 cells in thalamus
lowest_count = 10
mask = np.array([xx[0] for xx in cell_counts_per_brain_left[:, :1] > lowest_count])

_dcontra_hc = _dcontra.T[mask].T #remove thalamus soi
_dipsi_hc = _dipsi.T[mask].T
_ccontra_hc = _ccontra.T[mask].T
_cipsi_hc = _cipsi.T[mask].T
_ratio_hc = _ratio.T[mask].T

filter_dcontra = np.array([struct[_ratio_hc[0]>threshold] for struct in _dcontra_hc])
filter_dipsi = np.array([struct[_ratio_hc[0]>threshold] for struct in _dipsi_hc])
filter_ccontra = np.array([struct[_ratio_hc[0]>threshold] for struct in _ccontra_hc])
filter_cipsi = np.array([struct[_ratio_hc[0]>threshold] for struct in _cipsi_hc])

filter_primary_pool = primary_pool_filtered[_ratio_hc[0]>threshold]
filter_ak_pool = ak_pool[np.unique(filter_primary_pool)]
filter_primary_lob_n = np.asarray([np.where(filter_primary_pool == i)[0].shape[0] 
                                    for i in np.unique(filter_primary_pool)])
#show contra, ipsi, and contra+ipsi heatmaps side by side

#mean percent counts
_pccontra = np.nan_to_num(np.array([xx[1:]/xx[0] for xx in filter_ccontra.T])*100) #note that the whole thalamus is the first element in nthe array
mean_contra_pcounts = np.array([np.median(_pccontra[np.where(filter_primary_pool == idx)[0]], axis=0) 
                for idx in np.unique(filter_primary_pool)])

#get acronyms of nuclei
ann_df = pd.read_excel(df_pth)
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
#choose whether to annotate the numbers in the heatmap
annotate = False
#set label coords
xaxis_label_x,xaxis_label_y = 0.7, 0.1
#set range of colormap
vmin = 0
vmax = 7

fig, axes = plt.subplots(ncols = 3, nrows = 1, figsize = (8,6), sharey = True, gridspec_kw = {"wspace":0, "hspace":0})

ax = axes[0]
show = mean_contra_pcounts.T 

cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,((vmax-vmin)/2)+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.1, aspect=10)
cb.set_label("% of thalamic cells", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
if annotate:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col < 1:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        

# xticks
ax.set_xticks(np.arange(len(filter_ak_pool))+.5)
lbls = np.asarray(filter_ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, filter_primary_lob_n)], 
                    rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="small")
#make x label
ax.set_xlabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(xaxis_label_x,xaxis_label_y)

#ipsi side
ax = axes[1]

_pcipsi = np.array([xx[1:]/xx[0] for xx in filter_cipsi.T])*100
mean_ipsi_pcounts = np.asarray([np.median(_pcipsi[np.where(filter_primary_pool == idx)[0]], axis=0) 
            for idx in np.unique(filter_primary_pool)])

show = mean_ipsi_pcounts.T 

#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.1, aspect=10)
cb.set_label("% of thalamic cells", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
if annotate:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col < 1:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
# xticks
ax.set_xticks(np.arange(len(filter_ak_pool))+.5)
lbls = np.asarray(filter_ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, filter_primary_lob_n)], 
                    rotation=30, fontsize="x-small", ha="right")
#make x label
ax.set_xlabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(xaxis_label_x,xaxis_label_y)

#combined sides
#ipsi side
ax = axes[2]

_pcounts = np.array([xx[1:]/xx[0] for xx in (filter_cipsi+filter_ccontra).T])*100
mean_pcounts = np.asarray([np.mean(_pcounts[np.where(filter_primary_pool == idx)[0]], axis=0) for idx 
                           in np.unique(filter_primary_pool)])

show = mean_pcounts.T 

#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.1, aspect=10)
cb.set_label("% of thalamic cells", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
if annotate:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col < 1:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
# xticks
ax.set_xticks(np.arange(len(filter_ak_pool))+.5)
lbls = np.asarray(filter_ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, filter_primary_lob_n)], 
                    rotation=30, fontsize="x-small", ha="right")
#make x label
ax.set_xlabel("Bilateral", fontsize="small")
ax.yaxis.set_label_coords(xaxis_label_x,xaxis_label_y)

plt.savefig(os.path.join(fig_dst, "thal_contra_n_ipsi_mean_pcounts_threshold_%s.jpg" % threshold), 
            bbox_inches = "tight", dpi = 300)

plt.close()

#%%
#filter results by contra/ipsi ratio > 1 in thalamus

#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
pcounts_descending_order = np.sort(_pccontra)
order = np.argsort(np.median(_pccontra, axis = 0))
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
sois_descending_pcounts = np.array(short_nuclei)[order]

#boxplots of percent counts - contra
fig, axes = plt.subplots(ncols = 3, nrows = 1, figsize = (15,8), sharey = False, sharex = True, gridspec_kw = {"wspace":0, "hspace":0})

#first boxplot
ax = axes[0]

ax.boxplot(pcounts_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(pcounts_descending_order.T)
for i in range(ngroup):
    ax.scatter(pcounts_descending_order[:,i], 
                y=np.ones(len(pcounts_descending_order[:,0]))*i+1, color = "k", s = 10)


ax.set_xlabel("% of thalamic cells\nContra")
ax.set_ylabel("Thalamic nuclei")

#second boxplot, ipsi
ax = axes[1]

pcounts_descending_order = np.sort(_pcipsi)
order = np.argsort(np.median(_pcipsi, axis = 0))
sois_descending_pcounts = np.array(short_nuclei)[order]

ax.boxplot(pcounts_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(pcounts_descending_order.T)
for i in range(ngroup):
    ax.scatter(pcounts_descending_order[:,i], 
                y=np.ones(len(pcounts_descending_order[:,0]))*i+1, color = "k", s = 10)


ax.set_xlabel("Ipsi")

#final boxplot, combined
ax = axes[2]

pcounts_descending_order = np.sort(_pcounts)
order = np.argsort(np.median(_pcounts, axis = 0))
sois_descending_pcounts = np.array(short_nuclei)[order]

ax.boxplot(pcounts_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(pcounts_descending_order.T)
for i in range(ngroup):
    ax.scatter(pcounts_descending_order[:,i], 
                y=np.ones(len(pcounts_descending_order[:,0]))*i+1, color = "k", s = 10)

#label which boxplot belongs to which side
ax.set_xlabel("Bilateral")

plt.savefig(os.path.join(fig_dst, "thal_contra_n_ipsi_pcounts_boxplots_threshold_%s.pdf" % threshold), bbox_inches = "tight")

plt.close()
#%%
#do a rank correlation between nuclei for contra vs. ipsi side

ak_vh = np.array(["Vermis", "Hemisphere"])
func = lambda xx: 0 if xx < 2 else 1
filter_primary_pool_vh = np.array([func(xx) for xx in filter_primary_pool])

_pccontra_vermis = _pccontra[np.where(filter_primary_pool_vh == 0)]
_pccontra_hem = _pccontra[np.where(filter_primary_pool_vh == 1)]
_pcipsi_vermis = _pcipsi[np.where(filter_primary_pool_vh == 0)]
_pcipsi_hem = _pcipsi[np.where(filter_primary_pool_vh == 1)]

#vermis
contra_pcounts_descending_order = np.sort(_pccontra_vermis)
contra_order = np.argsort(np.mean(_pccontra_vermis, axis = 0))
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
contra_sois_descending_pcounts = np.array(short_nuclei)[contra_order]
contra_sois_ascending_pcounts = contra_sois_descending_pcounts[::-1]

ipsi_pcounts_descending_order = np.sort(_pcipsi_vermis)
ipsi_order = np.argsort(np.mean(_pcipsi_vermis, axis = 0))
ipsi_sois_descending_pcounts = np.array(short_nuclei)[ipsi_order]
ipsi_sois_ascending_pcounts = ipsi_sois_descending_pcounts[::-1]

contra_ranks = [i for i,nuc in enumerate(contra_sois_ascending_pcounts)]
ipsi_ranks = [i for j,nu in enumerate(ipsi_sois_ascending_pcounts) 
    for i,nuc in enumerate(contra_sois_ascending_pcounts) if ipsi_sois_ascending_pcounts[j] == nuc]

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#size of scatter
size = 70

Y = contra_ranks
X = ipsi_ranks

results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
ax.scatter(y = Y, x = X, s = size)

#plot fit line
ax.plot(mean_slope*range(len(X))+mean_intercept, '--k')    
    
ax.set_xlabel("Contralateral thalamus rank order (vermis injections)")
ax.set_ylabel("Ipsilateral thalamus rank order (vermis injections)")

#make text box
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)


plt.savefig(os.path.join(fig_dst, "thal_contra_v_ipsi_rank_order_vermis_threshold_%s.pdf" % threshold), bbox_inches = "tight")

plt.close()

#hemisphere
contra_pcounts_descending_order = np.sort(_pccontra_hem)
contra_order = np.argsort(np.mean(_pccontra_hem, axis = 0))
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
contra_sois_descending_pcounts = np.array(short_nuclei)[contra_order]
contra_sois_ascending_pcounts = contra_sois_descending_pcounts[::-1]

ipsi_pcounts_descending_order = np.sort(_pcipsi_hem)
ipsi_order = np.argsort(np.mean(_pcipsi_hem, axis = 0))
ipsi_sois_descending_pcounts = np.array(short_nuclei)[ipsi_order]
ipsi_sois_ascending_pcounts = ipsi_sois_descending_pcounts[::-1]

contra_ranks = [i for i,nuc in enumerate(contra_sois_ascending_pcounts)]
ipsi_ranks = [i for j,nu in enumerate(ipsi_sois_ascending_pcounts) 
    for i,nuc in enumerate(contra_sois_ascending_pcounts) if ipsi_sois_ascending_pcounts[j] == nuc]

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#size of scatter
size = 70

Y = contra_ranks
X = ipsi_ranks

results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
ax.scatter(y = Y, x = X, s = size)

#plot fit line
ax.plot(mean_slope*range(len(X))+mean_intercept, '--k')    
    
ax.set_xlabel("Contralateral thalamus rank order (hemisphere injections)")
ax.set_ylabel("Ipsilateral thalamus rank order (hemisphere injections)")

#make text box
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)


plt.savefig(os.path.join(fig_dst, "thal_contra_v_ipsi_rank_order_hemisphere_threshold_%s.pdf" % threshold), bbox_inches = "tight")

plt.close()

#%%
fig, axes = plt.subplots(ncols = 1, nrows = 6, figsize = (18,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [3,1.6,0.8,0.8,0.8,0.4]})


#set colormap specs
vmaxcounts = 90
whitetext = 10

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
ax.set_yticks(np.arange(len(ak))+.5)
ax.set_yticklabels(np.flipud(ak), fontsize="x-small")


ax = axes[1]
show = np.array([(density_per_brain_left[:, 2]+density_per_brain_right[:, 13])[np.argsort(_dist, axis = 0)], #vpl
                 (density_per_brain_left[:, 3]+density_per_brain_right[:, 13])[np.argsort(_dist, axis = 0)], #vpm
                 (density_per_brain_left[:, 10]+density_per_brain_right[:, 13])[np.argsort(_dist, axis = 0)], #md
                 (density_per_brain_left[:, 13]+density_per_brain_right[:, 13])[np.argsort(_dist, axis = 0)]]) #rtn
yaxis = ["VPL density", "VPM density", "MD density", "RTN density"]

vmin = 0
vmax = 20
cmap = plt.cm.Greens
cmap.set_over("darkgreen")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.8, aspect=10)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col > 15:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")

ax = axes[2]
show = sort_ccontra_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcounts
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.8, aspect=10)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)


ax = axes[3]
show = sort_cipsi_pool.T
yaxis = grps

vmin = 0
vmax = vmaxcounts
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.8, aspect=10)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)


ax = axes[4]
show = sort_dratio_pool.T
yaxis = grps

vmin = 0.7
vmax = 1.5
cmap = plt.cm.Blues
cmap.set_over("navy")
#colormap
bounds = np.linspace(vmin,vmax,5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
cb.set_label("Ratio", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col > 1.2:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra/Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[5]
show = np.asarray([sort_dist])
br = lr_brains 

vmin = -100
vmax = 80
cmap = plt.cm.RdBu_r
cmap.set_over('maroon')
cmap.set_under('midnightblue')
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=2, aspect=10)
cb.set_label("Left to right", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < -75 or col > 70:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")        

#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(sort_brains))+.5)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=5, ha="right")

ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

plt.savefig(os.path.join(fig_dst, "rtn_counts_ratios.pdf"), bbox_inches = "tight")
#%%
#basic statistics for these ratios

df = pd.DataFrame()
sort_dratio_pool[sort_dratio_pool == np.inf] = np.nan
mask_brains = np.ones(len(sort_dratio_pool)).astype(bool)
np.put(mask_brains, [1,3,6,7,16], np.zeros(len([1,3,6,7,16])).astype(bool))

from sklearn.impute import SimpleImputer

rat = SimpleImputer().fit_transform(sort_dratio_pool[mask_brains])

df["median ratio"] = np.median(rat, axis = 0)
df["mean ratio"] = np.mean(rat, axis = 0)
df["std ratio"] = np.std(rat, axis = 0)
df["est std ratio"] = mad(rat, axis = 0)/0.6745

df.index = grps
df = df.round(2)
df.to_csv(os.path.join(fig_dst, "thal_ratio_stats.csv"))