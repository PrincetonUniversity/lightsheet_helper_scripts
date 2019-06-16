#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:18:29 2019

@author: wanglab
"""

import os, pickle as pckl
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

data["mat"] = mat
data["pmat"] = pmat
data["mat_shuf"] = mat_shuf
data["p_shuf"] = p_shuf
data["ak"] = ak
data["primary_lob_n"] = primary_lob_n
data["regions"] = np.asarray(nuclei)
data["primary_pool"] = primary_pool

#store data (serialize)
with open(os.path.join(dst, "model_data.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#%%
#init dict
data = {}

data["brainnames"] = brains
data["expr_all_as_frac_of_lob"] = expr_all_as_frac_of_lob
data["expr_all_as_frac_of_inj"] = expr_all_as_frac_of_inj
data["primary_as_frac_of_lob"] = primary_as_frac_of_lob
data["secondary"] = secondary
data["cell_counts_per_brain"] = cell_counts_per_brain_pool
data["nc_regions"] = np.asarray(regions)
data["expr_all_as_frac_of_lob_pool"] = expr_all_as_frac_of_lob_pool
data["expr_all_as_frac_of_lob_pool_norm"] = normalised
data["expr_all_as_frac_of_inj_pool"] = expr_all_as_frac_of_inj_pool
data["primary_pool"] = primary
data["cb_regions_pool"] = ak_pool
data["cb_regions"] = ak

#store data (serialize)
with open(os.path.join(dst,"data.p"), "wb") as handle:
    pckl.dump(data, handle, protocol=pckl.HIGHEST_PROTOCOL)
    
#for json we have to de-pythonify a few things
#init dict
data = {}

data["brainnames"] = brains
data["expr_all_as_frac_of_lob"] = expr_all_as_frac_of_lob.tolist()
data["expr_all_as_frac_of_inj"] = expr_all_as_frac_of_inj.tolist()
data["secondary"] = secondary.tolist()
data["cell_counts_per_brain"] = cell_counts_per_brain_pool.tolist()
data["nc_regions"] = regions.tolist()
data["expr_all_as_frac_of_lob_pool"] = expr_all_as_frac_of_lob_pool.tolist()
data["expr_all_as_frac_of_lob_pool_norm"] = normalised.tolist()
data["expr_all_as_frac_of_inj_pool"] = expr_all_as_frac_of_inj_pool.tolist()
data["primary_pool"] = primary.tolist()
data["cb_regions_pool"] = ak_pool.tolist()
data["cb_regions"] = ak.tolist()

import json

#make it work for Python 2+3 and with Unicode
import io
try:
    to_unicode = unicode
except NameError:
    to_unicode = str

#write JSON file
with io.open(os.path.join(dst,"data.json"), 'w', encoding='utf8') as outfile:
    str_ = json.dumps(data,
                      indent=4, sort_keys=True,
                      separators=(',', ': '), ensure_ascii=False)
    outfile.write(to_unicode(str_))