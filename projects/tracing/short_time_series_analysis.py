#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:10:33 2020

@author: wanglab
"""

import os,pandas as pd,numpy as np,itertools

#dir
pth = "/jukebox/wang/zahra/tracing_projects/mapping_paper/revision_images/short_timepoint_counts/done"
#atlas ids
dcn_iids = [989,91,846] #these are the actual dcn ids, but i"m going to consider labels in the surrounding ids
#as dcn signal too (since i moved the atlas for better reg and the movement is not taken into account w the mapping)
#same principle for the other regions
vestnuc_iids = [209,202,225,217]
pons_iids = [931,574]
# dorsalcol_iids = [711,1039,903]
cuneate=711
gracile=1039
excuneate=903
#init empty dict to save the counts
dcn={}
dorcol = {}
cun={};gr={};excun={}
vn = {}
pons = {}
lst = os.listdir(pth); lst.sort()
for dfpth in lst:
    df = pd.read_csv(os.path.join(pth,dfpth))
    if "_1.csv" in dfpth:
        dcn[dfpth[:-6]] = len(df["Coordinate 1"].values)
    if "_2.csv" in dfpth:
        # dorcol[dfpth[:-6]] = len(df["Coordinate 1"].values)
        cun[dfpth[:-6]] = len(df[df["Segment IDs"]==cuneate]["Coordinate 1"])
        #add any additional cells to gracile
        gr[dfpth[:-6]] = len(df[df["Segment IDs"]==gracile]["Coordinate 1"])+len(df[~df["Segment IDs"].isin([cuneate,gracile,excuneate])]["Coordinate 1"])
        excun[dfpth[:-6]] = len(df[df["Segment IDs"]==excuneate]["Coordinate 1"])
    if "_3.csv" in dfpth:
        vn[dfpth[:-6]] = len(df["Coordinate 1"].values)
    if "_4.csv" in dfpth:
        pons[dfpth[:-6]] = len(df["Coordinate 1"].values)
    if "tp_bl6_ts04" in dfpth:
        #if it"s the old brain with nice reg
        dcn[dfpth[:-6]] = len(df[df["Segment IDs"].isin(dcn_iids)]["Coordinate 1"])
        # dorcol[dfpth[:-6]] = len(df[df["Segment IDs"].isin(dorsalcol_iids)]["Coordinate 1"])
        cun[dfpth[:-6]] = len(df[df["Segment IDs"]==cuneate]["Coordinate 1"])
        gr[dfpth[:-6]] = len(df[df["Segment IDs"]==gracile]["Coordinate 1"])
        excun[dfpth[:-6]] = len(df[df["Segment IDs"]==excuneate]["Coordinate 1"])    
        vn[dfpth[:-6]] = len(df[df["Segment IDs"].isin(vestnuc_iids)]["Coordinate 1"])
        pons[dfpth[:-6]] = len(df[df["Segment IDs"].isin(pons_iids)]["Coordinate 1"])

anndfpth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
anndf = pd.read_excel(anndfpth)
iids = dcn_iids+vestnuc_iids+pons_iids+[cuneate,gracile,excuneate]
vols = [anndf.loc[anndf["id"] == iid, "voxels_in_structure"].values[0] for iid in iids]
vols = [sum(vols[:3]),sum(vols[3:7]),sum(vols[7:9]),vols[9],vols[10],vols[11]] #order is dcn,vn,pons,cun,gr,excun
scale = 0.020 #mm
dcn_d={}; cun_d={}; gr_d={}; excun_d={}; vn_d={}; pons_d={}
#convert to densities
for k,v in dcn.items():
    dcn_d[k] = v/(vols[0]*(scale**3))
for k,v in vn.items():
    vn_d[k] = v/(vols[1]*(scale**3))
for k,v in pons.items():
    pons_d[k] = v/(vols[2]*(scale**3))
for k,v in cun.items():
    cun_d[k] = v/(vols[3]*(scale**3))
for k,v in gr.items():
    gr_d[k] = v/(vols[4]*(scale**3))
for k,v in excun.items():
    excun_d[k] = v/(vols[5]*(scale**3))                

#export
df=pd.DataFrame.from_dict(dcn,orient="index")
df["Deep cerebellar nuclei"] = df[0]
df["Vestibular nuclei"] = list(vn.values())
df["BPN+NRTP"] = list(pons.values())
df["External cuneate n."] = list(excun.values())
df["Cuneate n."] = list(cun.values())
df["Gracile n."] = list(gr.values())
df=df.drop(columns=[0])
# df.to_csv(os.path.join(os.path.dirname(pth),"short_tp_counts.csv"))

#concantenate dorsal column nuclei into one column
df=pd.DataFrame()
df["counts"]=list(itertools.chain.from_iterable([dcn.values(),pons.values(),cun.values(),gr.values(),excun.values()]))
df["density"]=list(itertools.chain.from_iterable([dcn_d.values(),pons_d.values(),cun_d.values(),gr_d.values(),excun_d.values()]))
df["brain"]=list(itertools.chain.from_iterable([dcn.keys(),
             pons.keys(),cun.keys(),gr.keys(),excun.keys()]))
df["structure"]=list(itertools.chain.from_iterable([["cerebellar nuclei"]*len(dcn.values()),
                ["pons"]*len(dcn.values()),
                ["cuneate"]*len(gr.values()),
                ["gracile"]*len(cun.values()),["external cuneate"]*len(excun.values())]))
#ratios
lst=["cuneate","gracile","external cuneate"]
#only for hsv
dfhsv=df[df.brain.str.match("HSV")]
dfprv=df[df.brain.str.match("PRV")]
hsvdcn_to_dorcol=np.median(dfhsv[dfhsv.structure.isin(lst)]["density"].values)/np.median(dfhsv[dfhsv.structure=="cerebellar nuclei"]["density"].values)
#sem
from scipy.stats import median_abs_deviation as mad
sigma=mad(dfhsv[dfhsv.structure.isin(lst)]["density"].values)
n=5#degrees of freedom same as # brains
dorcol_sigma_med=(sigma/0.6745)/np.sqrt(n)
sigma=mad(dfhsv[dfhsv.structure=="cerebellar nuclei"]["density"].values)
n=5
dcn_sigma_med=(sigma/0.6745)/np.sqrt(n)

#ratio
sem_hsvdcn_to_dorcol=dorcol_sigma_med/dcn_sigma_med
#%%

df=pd.DataFrame.from_dict(dcn_d,orient="index")
df["Deep cerebellar nuclei"] = df[0]
df["Vestibular nuclei"] = list(vn_d.values())
df["BPN+NRTP"] = list(pons_d.values())
df["External cuneate n."] = list(excun_d.values())
df["Cuneate n."] = list(cun_d.values())
df["Gracile n."] = list(gr_d.values())
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