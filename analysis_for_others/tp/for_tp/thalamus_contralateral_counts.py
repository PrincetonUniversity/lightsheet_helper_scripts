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

lr_brains = list(lr_dist.keys())
id_table = pd.read_excel(df_pth)

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


#%%
dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"

contra_df = pd.DataFrame(contra)
contra_df.to_csv(os.path.join(dst, "thal_contra_counts_23_brains.csv")) 

ipsi_df = pd.DataFrame(contra)
ipsi_df.to_csv(os.path.join(dst, "thal_ipsi_counts_23_brains.csv")) 

data["contra"] = contra
data["ipsi"] = ipsi
data["lr_brains"] = lr_brains
#store data (serialize)
with open(os.path.join(dst, "thal_contra_ipsi_counts_23_brains.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)


#############################################################################################################################################

#%%
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
