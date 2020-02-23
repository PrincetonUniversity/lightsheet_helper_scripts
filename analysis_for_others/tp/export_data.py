#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:18:29 2019

@author: wanglab
"""

import os, pickle as pckl

#%%


dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/"

data = {}

data["brains"] = brains
data["ak_pool"] = ak_pool
data["primary_lob_n"] = primary_lob_n
data["expr_all_as_frac_of_lob"] = expr_all_as_frac_of_lob
data["primary_as_frac_of_lob"] = primary_as_frac_of_lob
data["expr_all_as_frac_of_inj"] = expr_all_as_frac_of_inj

#store data (serialize)
with open(os.path.join(dst, "nc_hsv_maps_contra_pma.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)

#%%

#save
dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/"
data = {}

data["brainnames"] = brains
data["frac_of_inj_pool"] = frac_of_inj_pool
data["primary_pool"] = primary_pool
data["ak_pool"] = ak_pool
data["density_per_struct"] = density_per_struct
data["counts_per_struct"] = counts_per_struct
data["short_nuclei"] = short_nuclei
data["primary_lob_n"] = primary_lob_n
#store data (serialize)
with open(os.path.join(dst, "vtasnc_counts_contra.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)


#%%

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/"

data = {}

data["c_mat"] = c_mat
data["mat"] = mat
data["pmat"] = pmat
data["mat_shuf"] = mat_shuf
data["p_shuf"] = p_shuf
data["ak_pool"] = ak_pool
data["primary_lob_n"] = primary_lob_n
data["regions"] = regions
data["primary_pool"] = primary_pool
data["pcounts"] = pcounts
data["expr_all_as_frac_of_inj_pool_norm"] = expr_all_as_frac_of_inj_pool_norm
#store data (serialize)
with open(os.path.join(dst, "thal_model_data_contra_allen.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#%%  

dst = "/home/wanglab/Desktop/zahra"#/jukebox/wang/zahra/h129_contra_vs_ipsi/data/"

data = {}

data["c_mat"] = c_mat
data["mat"] = mat
data["pmat"] = pmat
data["mat_shuf"] = mat_shuf
data["p_shuf"] = p_shuf
data["ak_pool"] = ak_pool
data["primary_lob_n"] = primary_lob_n
data["regions"] = regions
data["primary_pool"] = primary_pool
data["pcounts_pool"] = pcounts_pool
data["frac_of_inj_pool"] = frac_of_inj_pool
#store data (serialize)
with open(os.path.join(dst, "nc_model_data_contra_pma.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#%%

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/"

data = {}

data["brains"] = brains
data["frac_of_inj_pool"] = frac_of_inj_pool
data["primary_pool"] = primary_pool
data["ak_pool"] = ak_pool
data["primary_lob_n"] = primary_lob_n


#store data (serialize)
with open(os.path.join(dst, "thal_hsv_maps_contra_allen.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#%%  

dst = "/jukebox/wang/zahra/tracing_projects/prv/for_tp/"

data = {}

data["c_mat"] = c_mat
data["mat"] = mat
data["pmat"] = pmat
data["mat_shuf"] = mat_shuf
data["p_shuf"] = p_shuf
data["ak_pool"] = ak_pool
data["primary_lob_n"] = primary_lob_n
data["regions"] = regions
data["primary_pool"] = primary_pool
data["pcounts_pool"] = pcounts_pool
data["frac_of_inj_pool"] = frac_of_inj_pool
#store data (serialize)
with open(os.path.join(dst, "model_data_contra_pma.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#%%

dst = "/jukebox/wang/zahra/tracing_projects/prv/for_tp/"

data = {}

data["brains"] = brains
data["frac_of_inj_pool"] = frac_of_inj_pool
data["primary_pool"] = primary_pool
data["ak_pool"] = ak_pool

#store data (serialize)
with open(os.path.join(dst, "prv_maps_contra_pma.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)

    
#%%

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
data = {}

data["lr_dist"] = lr_dist
data["thal_inj_vol"] = thal_inj_vol
#store data (serialize)
with open(os.path.join(dst, "thal_contra_counts.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)


#%%

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
data = {}

data["lr_dist"] = lr_dist
data["thal_inj_vol"] = thal_inj_vol
#store data (serialize)
with open(os.path.join(dst, "thal_inj_contra_v_ipsi.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)

#%%

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
data = {}

data["nc_inj_vol"] = nc_inj_vol
#store data (serialize)
with open(os.path.join(dst, "nc_inj_vol.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)

data["thal_inj_vol"] = thal_inj_vol
#store data (serialize)
with open(os.path.join(dst, "thal_inj_vol.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)



#%%

dst = "/jukebox/wang/zahra/h129_qc/data"
data = {}

data["nc_brainnames"] = nc_brains
data["thal_brainnames"] = thal_brains
data["thal_cell_counts_per_brain"] = thal_cell_counts_per_brain
data["nc_cell_counts_per_brain"] = nc_cell_counts_per_brain
data["nc_volume_per_brain"] = nc_volume_per_brain
data["thal_volume_per_brain"] = thal_volume_per_brain
data["thal_density_per_brain"] = thal_density_per_brain
data["mean_thal_density_per_brain"] = mean_thal_density_per_brain
data["std_thal_density_per_brain"] = std_thal_density_per_brain
data["nc_density_per_brain"] = nc_density_per_brain
data["mean_nc_density_per_brain"] = mean_nc_density_per_brain
data["std_nc_density_per_brain"] = std_nc_density_per_brain
data["nc_lbls"] = lbls

#store data (serialize)
with open(os.path.join(dst, "nc_density_at_thal_nc_timepoint_data_all_brains.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)

#%%
#save
dst = "/jukebox/wang/zahra/modeling/h129/striatum/"

data = {}

data["brainnames"] = brains
data["expr_all_as_frac_of_lob"] = expr_all_as_frac_of_lob
data["expr_all_as_frac_of_inj"] = expr_all_as_frac_of_inj
data["primary_as_frac_of_lob"] = primary_as_frac_of_lob
data["secondary"] = secondary
data["cell_counts_per_brain"] = cell_counts_per_brain
data["cell_counts_per_brain_p"] = cell_counts_per_brain_p
data["expr_all_as_frac_of_lob_pool"] = expr_all_as_frac_of_lob_pool
#data["expr_all_as_frac_of_lob_pool_norm"] = expr_all_as_frac_of_lob_pool_norm
data["expr_all_as_frac_of_inj_pool"] = expr_all_as_frac_of_inj_pool
data["primary_pool"] = primary_pool
data["ak_pool"] = ak_pool
data["volume_per_brain"] = volume_per_brain
data["density_per_brain"] = density_per_brain
data["primary_lob_n"] = primary_lob_n
data["sois"] = sois

#store data (serialize)
with open(os.path.join(dst, "count_and_density_data.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#%%

#save
dst = "/jukebox/wang/zahra/modeling/h129/pallidum/"

data = {}

data["brainnames"] = brains
data["expr_all_as_frac_of_lob"] = expr_all_as_frac_of_lob
data["expr_all_as_frac_of_inj"] = expr_all_as_frac_of_inj
data["primary_as_frac_of_lob"] = primary_as_frac_of_lob
data["secondary"] = secondary
data["cell_counts_per_brain"] = cell_counts_per_brain
data["cell_counts_per_brain_p"] = cell_counts_per_brain_p
data["expr_all_as_frac_of_lob_pool"] = expr_all_as_frac_of_lob_pool
#data["expr_all_as_frac_of_lob_pool_norm"] = expr_all_as_frac_of_lob_pool_norm
data["expr_all_as_frac_of_inj_pool"] = expr_all_as_frac_of_inj_pool
data["primary_pool"] = primary_pool
data["ak_pool"] = ak_pool
data["volume_per_brain"] = volume_per_brain
data["density_per_brain"] = density_per_brain
data["primary_lob_n"] = primary_lob_n
data["sois"] = sois

#store data (serialize)
with open(os.path.join(dst, "count_and_density_data.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
    
#%%
#save
dst = "/jukebox/wang/zahra/modeling/h129/thalamus/"

data = {}

data["cell_counts_per_brain_p"] = cell_counts_per_brain_p
data["fit"] = fit
data["fit_shuf"] = fit_shuf
data["p_shuf"] = p_shuf
data["brains"] = brains
data["ak_pool"] = ak_pool
data["primary_pool"] = primary_pool

#store data (serialize)
with open(os.path.join(dst, "shuffle_figure_data.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#%%
#save
dst = "/jukebox/wang/zahra/modeling/h129/neocortex/"

data = {}

data["cell_counts_per_brain_p"] = cell_counts_per_brain_pool
data["fit"] = fit
data["fit_shuf"] = fit_shuf
data["p_shuf"] = p_shuf
data["brains"] = brains
data["ak_pool"] = ak
data["primary_pool"] = primary_pool
data["regions"] = regions

#store data (serialize)
with open(os.path.join(dst, "shuffle_figure_data.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)

    