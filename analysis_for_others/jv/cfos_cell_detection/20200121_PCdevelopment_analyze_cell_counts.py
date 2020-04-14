#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 18:17:23 2019

@author: wanglab
"""


import os, numpy as np, pandas as pd, scipy, itertools, sys
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.analysis.network_analysis import make_structure_objects

src = "/jukebox/wang/Jess/lightsheet_output/201812_PCdevelopment/pooled_analysis/all_cfos"
flds = "/jukebox/wang/Jess/lightsheet_output/201812_PCdevelopment/forebrain/processed"

#get files
lst = [os.path.join(flds, fld) for fld in os.listdir(flds) if os.path.exists(os.path.join(os.path.join(flds, fld), 
                                                     "Annotated_counts_60um_edge_80um_vntric_erosion.csv"))]; lst.sort()
#conditions
nms = ["an_01_lob6", "an_02_lob6", "an_03_lob6", "an_04_lob6",
       "an_05_lob6", "an_06_lob6", "an_07_lob6", "an_09_lob6",
       "an_13_noinj", "an_14_noinj", "an_15_noinj", "an_16_noinj",
       "an_17_crus1", "an_18_crus1", "an_19_crus1", "an_20_crus1",
       "an_21_crus1", "an_22_crus1", "an_23_crus1", "an_24_crus1",
       "an_25_crus1", "an_26_crus1", "an_27_crus1"]

cond = ["JuvenilepretreatCNOandlobuleVI_VecCtrl", "JuvenilepretreatCNOandlobuleVI_DREADDs","JuvenilepretreatCNOandlobuleVI_VecCtrl",
        "JuvenilepretreatCNOandlobuleVI_DREADDs", "JuvenilepretreatCNOandlobuleVI_VecCtrl","JuvenilepretreatCNOandlobuleVI_DREADDs",
        "JuvenilepretreatCNOandlobuleVI_DREADDs", "JuvenilepretreatCNOandlobuleVI_VecCtrl",
        "untreated","untreated","untreated","untreated", "JuvenilepretreatCNOandcrusI_DREADDs","JuvenilepretreatCNOandcrusI_VecCtrl",
        "JuvenilepretreatCNOandcrusI_VecCtrl","JuvenilepretreatCNOandcrusI_DREADDs", "JuvenilepretreatCNOandcrusI_VecCtrl",
        "JuvenilepretreatCNOandcrusI_DREADDs", "JuvenilepretreatCNOandcrusI_DREADDs",
        "JuvenilepretreatCNOandcrusI_DREADDs","JuvenilepretreatCNOandcrusI_VecCtrl",
        "JuvenilepretreatCNOandcrusI_DREADDs","JuvenilepretreatCNOandcrusI_DREADDs"]
        
conditions = {n:c for n,c in zip(nms, cond)}

#set appropriate paths
pth = os.path.join(src, "cell_counts.csv")
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
curated_structures = "/jukebox/wang/Jess/lightsheet_output/201812_PCdevelopment/structures.csv"

#build structures class
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)


def generate_data_frame(conditions, lst, pth, flnm):
    """ 
    used to make a pooled csv file of all cell counts in an experiment
    inputs:
        conditions: zip of file names + condition
        lst: list of file names run through analysis
        pth: path to save csv output
    """
    #generate data frame
    bglst=[]
    for fl in lst:
        #extract out info
        nm = os.path.basename(fl)
        #make dataframe
        df = pd.read_csv(fl+"/%s.csv" % flnm)[1:] #remove previous headers
        print(nm, df.shape)
        df = df.replace(np.nan, "", regex=True)
        df["Brain"] = nm
        df["Condition"] = conditions[nm]
        bglst.append(df)
    
    df = pd.concat(bglst)
    df["counts"] = df["counts"].apply(int)

    #remove structures we don"t care about!!!!!!!!!!!!!!!!!!!!!!
    sois = ["ventricular systems", "fiber tracts", "grooves"]
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        df = df[df.name != soi.name]
        progeny = [str(xx.name) for xx in soi.progeny]
        for progen in progeny:
            df = df[df.name != progen]
    
    df.drop(columns = ["Unnamed: 0"]).to_csv(pth, index = None)
    
    return pth

#run
csv_pth = generate_data_frame(conditions, lst, pth, "Annotated_counts_60um_edge_80um_vntric_erosion")


def generate_percent_counts_and_density_per_region(src, csv_pth):
    
    """ generates another column in the dataframe that is just # counts / total counts in brain and density per region"""
    
    #read csv
    df = pd.read_csv(csv_pth)
    
    #set resolution of atlas
    scale = 0.020 ##mm/voxel
    
    #get all brain names
    brains = np.unique(df.Brain.values)
    
    #save total counts in dict
    total_counts = {}
    percents = {}
    #for each brain, get total counts
    for brain in brains:
        total_counts[brain] = df[df.Brain == brain].counts.sum(0)    
           
    percents = [df[df.Brain == brain].counts.apply(lambda x: (x/total_counts[brain])*100).astype("float64") for brain in brains]
    
    #concantenate together
    df_percent = pd.concat(percents)
    
    df["percent"] = df_percent
    df["density"] = df[df.voxels_in_structure > 0].apply(lambda x:x.counts/(float(x.voxels_in_structure*(scale**3))), 1)
    
    #export
    df.to_csv(os.path.join(src, "cell_counts.csv"))
    
    return os.path.join(src, "cell_counts.csv")

#run
percent_density_csv_pth = generate_percent_counts_and_density_per_region(src, csv_pth)

#do stats for previously curated structures
sois = pd.read_csv(curated_structures)
sois = [xx[0] for xx in sois.values]

def pool_regions(sois, percent_density_csv_pth, measures = ["counts", "percent", "density"], scale = 0.020):
    #set variables
    orgdf = pd.read_csv(percent_density_csv_pth, index_col = None)
    #init
    tdf = orgdf.copy()
    #we have to do this per animal
    #get the brains
    brains = tdf.Brain.unique()
    #init giant df
    main_df = pd.DataFrame()

    for an in brains:
        for m, measure in enumerate(measures): #iterate through the measures
            df = tdf[tdf.Brain == an]
            condition = df["Condition"].unique()
            pooled_counts = {}
            
            #loop through and find downstream structures of a given parent
            for soi in sois:
                soi = [s for s in structures if s.name==soi][0]
                progeny = [str(xx.name) for xx in soi.progeny]
                counts = [] #store counts in this list
                volume = []
                val = df.loc[(df.name == soi.name), measure].values
                vol = df.loc[(df.name == soi.name), "voxels_in_structure"].values
                if val.shape[0] > 0: counts.append(val[0])
                if vol.shape[0] > 0: volume.append(vol[0])
                if len(progeny) > 0:
                    for progen in progeny:
                        val = df.loc[(df.name == progen), measure].values
                        vol = df.loc[(df.name == progen), "voxels_in_structure"].values
                        if val.shape[0] > 0: counts.append(val[0])
                        if vol.shape[0] > 0: volume.append(vol[0])
                        
                #sum counts in parent structures            
                pooled_counts[soi.name] = np.sum(np.asarray(counts))
                if measure == "density": pooled_counts[soi.name] = np.sum(np.array(counts))/(float(np.sum(np.array(volume))*(scale**3)))
                
            #fill other details and add to big dict
            if m == 0:
                pooled_counts_an = pd.DataFrame(list(pooled_counts.values()))
                pooled_counts_an.columns = ["%s" % measure]
                pooled_counts_an["name"] = list(pooled_counts.keys())
            else:
                pooled_counts_an["%s" % measure] = list(pooled_counts.values())
            
        pooled_counts_an["Brain"] = list(itertools.repeat(an, len(pooled_counts_an)))
        pooled_counts_an["Condition"] = list(itertools.repeat(condition[0], len(pooled_counts_an)))
        main_df = main_df.append(pooled_counts_an)
        
    #save
    main_df.to_csv(os.path.join(src, "cell_counts_curated_structures.csv"))
    
    print("saved in :{}".format(src))
    
    return os.path.join(src, "cell_counts_curated_structures.csv")

curated_csv = pool_regions(sois, percent_density_csv_pth)

#t-stats
def generate_paired_statistics(src, csv_pth, conditions, measures = ["counts", "percent", "density"], curated = False):
    """
    generates paried t tests, and Mann Whitney and Wilcox Rank test results from pooled csv counts
    returns:
        tdf_dct: data frame with comparison tests
        sigs: significant structures (p < 0.05)
    """
    
    df = pd.read_csv(csv_pth)
    tdf_dct={}
    
    stim, ctrl = conditions
    df = df[df["Condition"].isin([stim, ctrl])]                                                                                                             
    structure_list = df.name.unique()
    brains = df.Brain.unique()
    
    print("*************")
    print(stim, ctrl, "len of brains: {}".format(len(brains)))
    for measure in measures: #iterate thru measures
        lst = []
        
        for structure in structure_list:
            
                smean = np.float32(df.loc[((df.name == structure) & (df.Condition == stim)), measure].mean())
                sstd = np.float32(df.loc[((df.name == structure) & (df.Condition == stim)), measure].std())
                
                cmean = np.float32(df.loc[((df.name == structure) & (df.Condition == ctrl)), measure].mean())
                cstd = np.float32(df.loc[((df.name == structure) & (df.Condition == ctrl)), measure].std())
                pval = scipy.stats.ttest_ind(np.float32(df.loc[((df.name == structure) & (df.Condition == ctrl)), measure]), 
                                                     np.float32(df.loc[((df.name == structure) & (df.Condition == stim)), measure]))
                try:
                    mannwhit = scipy.stats.mannwhitneyu(np.float32(df.loc[((df.name == structure) & (df.Condition == ctrl)), measure]), 
                                                         np.float32(df.loc[((df.name == structure) & (df.Condition == stim)), measure]), alternative = "two-sided")
                except ValueError:
                    mannwhit = [0.0,1.0]
                
                wilcoxrank = scipy.stats.ranksums(np.float32(df.loc[((df.name == structure) & (df.Condition == ctrl)), measure]), 
                                                     np.float32(df.loc[((df.name == structure) & (df.Condition == stim)), measure]))
                lst.append((structure, smean, sstd, cmean, cstd, pval[0], pval[1], mannwhit[1], wilcoxrank[1]))
                
    
        #print lst
        #make tmp df
        tdf = pd.DataFrame(data=lst, columns=["name",  "%s mean" % stim, "%s std" % stim, 
                                              "%s mean" % ctrl, "%s std" % ctrl,"tstat", "pval", "mannwhit", "wilcoxrank"])
        tdf.sort_values("mannwhit")
        if curated:
            tdf.to_csv(os.path.join(src, "{}_mannwhit_{}-{}_curated.csv".format(measure, ctrl, stim)))
            print("saved to: {}".format(os.path.join(src, "{}_mannwhit_{}-{}_curated.csv".format(measure, ctrl, stim))))
        else:
            tdf.to_csv(os.path.join(src, "{}_mannwhit_{}-{}.csv".format(measure, ctrl, stim)))
            print("saved to: {}".format(os.path.join(src, "{}_mannwhit_{}-{}.csv".format(measure, ctrl, stim))))
        tdf_dct["{} vs {}".format(ctrl, stim)]=tdf #this vs is important down the road
            
    return

#run
tdf_dct = generate_paired_statistics(src, curated_csv, ["JuvenilepretreatCNOandlobuleVI_DREADDs", #format as [stimulation condition, control]
                                                        "JuvenilepretreatCNOandlobuleVI_VecCtrl"], ["percent"], curated = True)
    
tdf_dct = generate_paired_statistics(src, curated_csv, ["JuvenilepretreatCNOandcrusI_DREADDs", #format as [stimulation condition, control]
                                                        "JuvenilepretreatCNOandcrusI_VecCtrl"], ["percent"], curated = True)
    
tdf_dct = generate_paired_statistics(src, curated_csv, ["JuvenilepretreatCNOandcrusI_DREADDs", #format as [stimulation condition, control]
                                                        "untreated"], ["percent"], curated = True)

tdf_dct = generate_paired_statistics(src, curated_csv, ["JuvenilepretreatCNOandlobuleVI_DREADDs", #format as [stimulation condition, control]
                                                        "untreated"], ["percent"], curated = True)
        
#%%
#one way anova
from scipy.stats import f_oneway
#pooled results
df = pd.read_csv(curated_csv, index_col = None).drop(columns = ["Unnamed: 0"])

df_anova = pd.DataFrame()
df_anova["name"] = df.name.unique()

for nm in np.unique(df.name.values): #only gets unique names
    
    f, pval = f_oneway(df[(df.name == nm) & (df.Condition == "JuvenilepretreatCNOandlobuleVI_DREADDs")]["percent"].values, 
                 df[(df.name == nm) & (df.Condition == "JuvenilepretreatCNOandlobuleVI_VecCtrl")]["percent"].values,
                 df[(df.name == nm) & (df.Condition == "untreated")]["percent"].values)
    
    df_anova.loc[(df_anova["name"] == nm), "anova_percent_counts_f"] = f
    df_anova.loc[(df_anova["name"] == nm), "anova_percent_counts_pval"] = pval
        
df_anova.to_csv(os.path.join(src, "one_way_anova_lobVI_curated_structures.csv"))

#%%
#for crus I
df_anova = pd.DataFrame()
df_anova["name"] = df.name.unique()

for nm in np.unique(df.name.values): #only gets unique names
    
    f, pval = f_oneway(df[(df.name == nm) & (df.Condition == "JuvenilepretreatCNOandcrusI_DREADDs")]["percent"].values, 
                 df[(df.name == nm) & (df.Condition == "JuvenilepretreatCNOandcrusI_VecCtrl")]["percent"].values,
                 df[(df.name == nm) & (df.Condition == "untreated")]["percent"].values)
    
    df_anova.loc[(df_anova["name"] == nm), "anova_percent_counts_f"] = f
    df_anova.loc[(df_anova["name"] == nm), "anova_percent_counts_pval"] = pval
        
df_anova.to_csv(os.path.join(src, "one_way_anova_crusI_curated_structures.csv"))