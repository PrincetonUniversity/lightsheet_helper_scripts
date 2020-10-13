#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:10:33 2020

@author: wanglab
"""

import os,pandas as pd,numpy as np

#dir
pth = "/jukebox/wang/zahra/tracing_projects/mapping_paper/revision_images/short_timepoint_counts/done"
#atlas ids
dcn_iids = [989,91,846] #these are the actual dcn ids, but i"m going to consider labels in the surrounding ids
#as dcn signal too (since i moved the atlas for better reg and the movement is not taken into account w the mapping)
#same principle for the other regions
vestnuc_iids = [209,202,225,217]
pons_iids = [931,574]
dorsalcol_iids = [711,1039,903]
#init empty dict to save the counts
dcn = {}
dorcol = {}
vn = {}
pons = {}
lst = os.listdir(pth); lst.sort()
for dfpth in lst:
    df = pd.read_csv(os.path.join(pth,dfpth))
    if "_1.csv" in dfpth:
        dcn[dfpth[:-6]] = len(df["Coordinate 1"].values)
    if "_2.csv" in dfpth:
        dorcol[dfpth[:-6]] = len(df["Coordinate 1"].values)
    if "_3.csv" in dfpth:
        vn[dfpth[:-6]] = len(df["Coordinate 1"].values)
    if "_4.csv" in dfpth:
        pons[dfpth[:-6]] = len(df["Coordinate 1"].values)
    if "tp_bl6_ts04" in dfpth:
        #if it"s the old brain with nice reg
        dcn[dfpth[:-6]] = len(df[df["Segment IDs"].isin(dcn_iids)]["Coordinate 1"])
        dorcol[dfpth[:-6]] = len(df[df["Segment IDs"].isin(dorsalcol_iids)]["Coordinate 1"])
        vn[dfpth[:-6]] = len(df[df["Segment IDs"].isin(vestnuc_iids)]["Coordinate 1"])
        pons[dfpth[:-6]] = len(df[df["Segment IDs"].isin(pons_iids)]["Coordinate 1"])

anndfpth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
anndf = pd.read_excel(anndfpth)
iids = dcn_iids+vestnuc_iids+pons_iids+dorsalcol_iids
vols = [anndf.loc[anndf["id"] == iid, "voxels_in_structure"].values[0] for iid in iids]
vols = [sum(vols[:3]),sum(vols[3:7]),sum(vols[7:10]),sum(vols[10:])]
scale = 0.020 #mm
dcn_d={}; dorcol_d={}; vn_d={}; pons_d={}
#convert to densities
for k,v in dcn.items():
    dcn_d[k] = v/(vols[0]*(scale**3))
for k,v in dorcol.items():
    dorcol_d[k] = v/(vols[0]*(scale**3))
for k,v in vn.items():
    vn_d[k] = v/(vols[0]*(scale**3))
for k,v in pons.items():
    pons_d[k] = v/(vols[0]*(scale**3))

#export
df=pd.DataFrame.from_dict(dcn,orient="index")
df["Dorsal column nuclei"] = list(dorcol.values())
df["Vestibular nuclei"] = list(vn.values())
df["BPN+NRTP"] = list(pons.values())
df["Deep cerebellar nuclei"] = df[0]
df=df.drop(columns=[0])
df.to_csv(os.path.join(os.path.dirname(pth),"short_tp_counts.csv"))

df=pd.DataFrame.from_dict(dcn_d,orient="index")
df["Dorsal column nuclei"] = list(dorcol_d.values())
df["Vestibular nuclei"] = list(vn_d.values())
df["BPN+NRTP"] = list(pons_d.values())
df["Deep cerebellar nuclei"] = df[0]
df=df.drop(columns=[0])
df=df.round(2)
df.to_csv(os.path.join(os.path.dirname(pth),"short_tp_density_cellsmm3.csv"))

#calculate ratios of hsv
dfhsv = df[df.index.str.startswith("HSV")]
dfctb = df[df.index.str.startswith("CTB")]
dfhsv["ratio DCN/Dorsal column"] = dfhsv["Deep cerebellar nuclei"]/dfhsv["Dorsal column nuclei"]
dfhsv["ratio DCN/BPN+NRTP"] = dfhsv["Deep cerebellar nuclei"]/dfhsv["BPN+NRTP"]
dfhsv["ratio DCN/Vestibular nuclei"] = dfhsv["Deep cerebellar nuclei"]/dfhsv["Vestibular nuclei"]

# ratio_dcn_dorcol = np.array(list(dcn.values()))/np.array(list(dorcol.values()))
# ratio_dcn_pons = np.array(list(dcn.values()))/np.array(list(pons.values()))
# ratio_dcn_vn = np.array(list(dcn.values()))/np.array(list(vn.values()))
# print(ratio_dcn_dorcol.mean())
# print(ratio_dcn_pons.mean())
# print(ratio_dcn_vn.mean())