#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 17:34:06 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl
from scipy.stats import median_absolute_deviation as mad
from tools.utils.io import makedir
from tools.analysis.network_analysis import make_structure_objects
import matplotlib.colors

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop
 
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#set paths
dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"
fig_dst = "/home/wanglab/Desktop"

ann_pth = os.path.join(dst, "atlases/sagittal_allen_ann_25um_iso_60um_edge_80um_ventricular_erosion.tif")

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

#import data
data_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_inj_contra_v_ipsi.p"
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")
lr_dist = data["lr_dist"]
thal_inj_vol = data["thal_inj_vol"]

#make structures
df_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
ann_df = pd.read_excel(df_pth)

structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
lr_brains = list(lr_dist.keys())
atl_dst = os.path.join(dst, "pma_to_aba"); makedir(atl_dst)
id_table = pd.read_excel(df_pth)


def get_cell_n_density_counts(brains, structure, structures, cells_regions, scale_factor = 0.025):
    """ consolidating to one function bc then no need to copy/paste """
    #get cell counts adn densities
    #get densities for all the structures
    df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx", index_col = None)
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
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0]/2) #divide by 2 since these are half brains!!!
                progeny = [str(xx.name) for xx in soi.progeny]
                #now sum up progeny
                if len(progeny) > 0:
                    for progen in progeny:
                        for k, v in cells_regions[brain].items():
                            if k == progen and progen != "Primary somatosensory area, unassigned, layer 4,5,6":
                                counts.append(v)
                                #add to volume list from LUT
                                volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0]/2) #divide by 2 since these are half brains!!!
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
            except Exception as e:
                print(e)
                for k, v in cells_regions[brain].items():
                    if k == soi:
                        counts.append(v)                    
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0]/2) #divide by 2 since these are half brains!!!
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
cells_regions = pckl.load(open(os.path.join(dst, "thal_right_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

#nuclei = ["Pons", "Thalamus, sensory-motor cortex related", "Thalamus, polymodal association cortex related"]
nuclei = ["Thalamus", "Ventral posteromedial nucleus of the thalamus", "Ventral posterolateral nucleus of the thalamus",
          "Ventral anterior-lateral complex of the thalamus", "Ventral medial nucleus of the thalamus", "Anteroventral nucleus of thalamus", 
          "Reticular nucleus of the thalamus", "Ventral part of the lateral geniculate complex", "Mediodorsal nucleus of thalamus",
          "Submedial nucleus of the thalamus", "Nucleus of reuniens", "Paraventricular nucleus of the thalamus", 
          "Central lateral nucleus of the thalamus", "Parafascicular nucleus", "Posterior complex of the thalamus",
          "Lateral dorsal nucleus of thalamus", "Lateral posterior nucleus of the thalamus", "Lateral habenula"]
#nuclei = ["Ventral anterior-lateral complex of the thalamus", "Ventral medial nucleus of the thalamus", 
#          "Ventral posterolateral nucleus of the thalamus", "Ventral posteromedial nucleus of the thalamus",
#          "Subparafascicular nucleus", "Subparafascicular area",
#          "Peripeduncular nucleus", "Medial geniculate complex", "Dorsal part of the lateral geniculate complex",
#          "Lateral posterior nucleus of the thalamus", "Posterior complex of the thalamus",
#          "Posterior limiting nucleus of the thalamus", "Suprageniculate nucleus", "Anteroventral nucleus of thalamus", 
#          "Anteromedial nucleus", "Anterodorsal nucleus", "Interanteromedial nucleus of the thalamus", 
#          "Interanterodorsal nucleus of the thalamus", "Lateral dorsal nucleus of thalamus", 
#          "Intermediodorsal nucleus of the thalamus", "Mediodorsal nucleus of thalamus", "Submedial nucleus of the thalamus", 
#          "Perireunensis nucleus", "Paraventricular nucleus of the thalamus", "Parataenial nucleus", "Nucleus of reuniens", 
#          "Rhomboid nucleus", "Central medial nucleus of the thalamus", "Paracentral nucleus",
#          "Central lateral nucleus of the thalamus", "Parafascicular nucleus", 
#          "Reticular nucleus of the thalamus", "Ventral part of the lateral geniculate complex", "Epithalamus"]

#RIGHT SIDE
cell_counts_per_brain_right, density_per_brain_right, volume_per_brain_right = get_cell_n_density_counts(brains, nuclei, 
                                                                                                         structures, cells_regions)
#LEFT SIDE
cells_regions = pckl.load(open(os.path.join(dst, "thal_left_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      
cell_counts_per_brain_left, density_per_brain_left, volume_per_brain_left = get_cell_n_density_counts(brains, nuclei, 
                                                                                                      structures, cells_regions)

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
        
#############################################################################################################################################
#USE THESE ARRAY FOR THE CONTRA COUNTS + DENSITIES
#NOTE: orientation is nuclei (x) vs. brains (y)
_ccontra = np.asarray(_ccontra).T; _dcontra = np.asarray(_dcontra).T
_cipsi = np.asarray(_cipsi).T; _dipsi = np.asarray(_dipsi).T
#############################################################################################################################################

_ratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
#make into one
_dist = np.asarray(list(lr_dist.values()))

#injection site analysis
data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data.p"
model_data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/model_data_v2.p"
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")
model_data = pckl.load(open(model_data_pth, "rb"), encoding = "latin1")

_primary_pool = data["primary_pool"]
ak_pool = data["cb_regions_pool"]
_inj = data["expr_all_as_frac_of_inj_pool"]
primary_lob_n = model_data["primary_lob_n"]


#%%
#filter results by contra/ipsi ratio > 1 in thalamus
#show contra, ipsi, and contra+ipsi heatmaps side by side

#mean percent counts
_pccontra = np.nan_to_num(np.array([xx[1:]/xx[0] for xx in _ccontra.T])*100) #note that the whole thalamus is the first element in nthe array
mean_contra_pcounts = np.array([np.mean(_pccontra[np.where(_primary_pool == idx)[0]], axis=0) 
                        for idx in np.unique(_primary_pool)])

#get acronyms of nuclei
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
#choose whether to annotate the numbers in the heatmap
annotate = False
#set label coords
xaxis_label_x,xaxis_label_y = 0.7, 0.1
#set range of colormap
vmin = 0
vmax = 8

fig = plt.figure(figsize=(5,7))
ax = fig.add_axes([.4,.1,.5,.8])

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
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="small")
#make x label
ax.set_xlabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(xaxis_label_x,xaxis_label_y)

plt.savefig(os.path.join(fig_dst, "thal_contra_mean_pcounts.pdf"), bbox_inches = "tight")

plt.close()

#%%

#basic statistics for these ratios
#NOTE: this is only if you want to look at sm vs. poly thal and collect the cell counts as such
df = pd.DataFrame()
#decimal to round by
d = 2
df["median"] = np.round(np.median(_ratio, axis = 0), d)
df["mean"] = np.round(np.mean(_ratio, axis = 0), d)
df["std"] = np.round(np.std(_ratio, axis = 0), d)
df["est std"] = np.round(mad(_ratio, axis = 0)/0.6745, d)

df.index = ["Thalamus, sensory-motor cortex related", "Thalamus, polymodal association cortex related"]

df.to_csv(os.path.join(fig_dst, "thal_contra_ipsi_ratio_stats.csv"))
