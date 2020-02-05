#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:34:55 2020

@author: wanglab
"""

import pandas as pd, os
os.chdir("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.analysis.network_analysis import make_structure_objects

pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
dst = "/jukebox/wang/zahra/registration_error_pma"

#make structures to traverse hierarchy
structures = make_structure_objects(pth)

#read main annotation LUT
anndf = pd.read_excel(pth)

#make new df to only save child structures that don't belong to NC, white matter, or ventricles
sois = ["Isocortex", "ventricular systems", "fiber tracts", "grooves"]
for soi in sois:
    soi = [s for s in structures if s.name==soi][0]
    anndf = anndf[anndf.name != soi.name]
    progeny = [str(xx.name) for xx in soi.progeny]
    for progen in progeny:
        anndf = anndf[anndf.name != progen]
        
#format    
anndf = anndf.drop(columns = ["Unnamed: 0", "cell_count"])
anndf.to_csv(os.path.join(dst, "LUT_child_structures_wo_NC_ventr_fiber_tracts.csv"))
