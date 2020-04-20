#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 20:18:36 2020

@author: wanglab
"""

import os, numpy as np, pandas as pd, scipy, itertools, sys, json, seaborn as sns
import matplotlib.pyplot as plt, matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

def get_progeny(dic,parent_structure,progeny_list):
    
    """ 
    ---PURPOSE---
    Get a list of all progeny of a structure name.
    This is a recursive function which is why progeny_list is an
    argument and is not returned.
    ---INPUT---
    dic                  A dictionary representing the JSON file 
                         which contains the ontology of interest
    parent_structure     The structure
    progeny_list         The list to which this function will 
                         append the progeny structures. 
    """
    if 'msg' in list(dic.keys()): dic = dic['msg'][0]
    
    name = dic.get('name')
    children = dic.get('children')
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get('name')
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get('name')
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 

def generate_data_frame(conditions, lst, pth, flnm, 
                        exclude, ontology_file, scale = 0.025):
    
    """ 
    used to make a pooled csv file of all cell counts in an experiment
    inputs:
        conditions: zip of file names + condition
        lst: list of file names run through analysis
        pth: path to save csv output
        scale: atlas scale, default 25 um/voxel (0.025 mm/voxel)
    """
    # #open ontology file
    with open(ontology_file) as json_file:
        ontology_dict = json.load(json_file)

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
        
        #sum up ontology for each structure
        for struct in df.name.values:
            progeny = []; counts = []
            get_progeny(ontology_dict, struct, progeny)
            #add structure counts to itself
            counts.append(df[df.name == struct]["counts"].values[0])
            for progen in progeny:
                try:
                    counts.append(df[df.name == progen]["counts"].values[0])
                except:
                    counts.append(0)
            df.loc[df.name == struct, "counts"] = np.sum(np.array(counts))
        
        #remove structures you don't want to analyze
        for soi in exclude:
            df = df[df.name != soi]
            progeny = []; get_progeny(ontology_dict, soi, progeny)
            for progen in progeny:
                df = df[df.name != progen]
        
        #append to composite dataframe
        bglst.append(df)

    df = pd.concat(bglst)
    df["counts"] = df["counts"].apply(int)

    df = df.drop(columns = ["Unnamed: 0"])
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
    
    df.to_csv(pth, index = None)
    
    return pth

#%%
if __name__ == "__main__":
      
    src = "/jukebox/LightSheetData/falkner-mouse/scooter/clearmap_processed"
    dst = "/jukebox/LightSheetData/falkner-mouse/scooter/pooled_analysis"
    
    #get files
    lst = [os.path.join(src, fld) for fld in os.listdir(src) 
           if os.path.exists(os.path.join(os.path.join(src, fld), 
                                                         "Annotated_counts.csv"))]; lst.sort()
    #conditions
    nms = ["mfnp3", "fmnp4", "mfnp2", "mmnp4", "mmnp6", 
           "mmnp5", "fmnp6", "fmnp5"]
    
    cond = ["male-female", "female-male", "male-female", "male-male", "male-male",
            "male-male", "female-male", "female-male"]
    
    conditions = {n:c for n,c in zip(nms, cond)}
    
    #set appropriate paths
    pth = os.path.join(dst, "cell_counts.csv")
    df_pth = "/home/wanglab/LightSheetData/falkner-mouse/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
    ann_pth = "/home/wanglab/LightSheetData/falkner-mouse/allen_atlas/annotation_2017_25um_sagittal_forDVscans.nrrd"
    atl_pth = "/home/wanglab/LightSheetData/falkner-mouse/allen_atlas/average_template_25_sagittal_forDVscans.tif"
    ontology_file = "/jukebox/LightSheetData/falkner-mouse/allen_atlas/allen.json"

    #name structures you want to exclude from allen hierarchy
    exclude = ["ventricular systems", "fiber tracts", "grooves"]
    
    #run
    csv_pth = generate_data_frame(conditions, lst, pth, "Annotated_counts", exclude, ontology_file)
    #%%
    #one way anova
    from scipy.stats import f_oneway
    from statsmodels.stats.multicomp import MultiComparison
    #pooled results
    
    df = pd.read_csv(csv_pth, index_col = None)
    
    df_anova = pd.DataFrame()
    df_anova["name"] = df.name.unique()
    df_anova["parent_name"] = [df.loc[df.name == nm, "parent_name"].values[0]
                               for nm in df_anova["name"]]
    df_anova["tukeyhsd"] = np.ones(len(df_anova))*np.nan

    #measure to use for tests
    measure = "density"
    for nm in np.unique(df.name.values): #only gets unique names
        
        f, pval = f_oneway(df[(df.name == nm) & (df.Condition == "male-male")][measure].values, 
                     df[(df.name == nm) & (df.Condition == "female-male")][measure].values,
                     df[(df.name == nm) & (df.Condition == "male-female")][measure].values)
        
        df_anova.loc[(df_anova["name"] == nm), "anova_f"] = f
        df_anova.loc[(df_anova["name"] == nm), "anova_pval"] = pval
        
        #doing post hoc on significant structures
        if pval < 0.05:
            mc = MultiComparison(df[df.name == nm][measure].values, df[df.name == nm ].Condition.values)
            result = mc.tukeyhsd(alpha=0.05) #cutoff at 0.05
            
            df_anova.loc[(df_anova["name"] == nm), "tukeyhsd"] = np.min(result.pvalues)
                
    df_anova.to_csv(os.path.join(dst, "one_way_anova_{}.csv".format(measure)))
#%%
    #find structures with multiple comparisons corrections    
    sig_structs = df_anova[(df_anova.anova_pval < 0.05) & 
                            (df_anova.tukeyhsd < 0.05)].name.values
    #plot the ones with lower density separately
    sig_structures_ld = [xx for xx in sig_structs 
                         if df[df.name == xx][measure].mean() < 100]
    
    df_nm = df[(df.name.isin(sig_structs)) & (~df.name.isin(sig_structures_ld))]
    fig, ax = plt.subplots(figsize=(5,10))
    sns.stripplot(y = "name", x = measure, hue = "Condition", data = df_nm,
                orient = "h", size=7)
    # sns.boxplot(y = "name", x = "density", hue = "Condition", data = df_nm,
    #             orient = "h", showfliers=False, showcaps=False, 
    #         boxprops={'facecolor':'None'})
    plt.xlabel("Density (cells/mm$^3$)")
    plt.ylabel("Structure")
    plt.savefig("/home/wanglab/Desktop/boxplots_%s.pdf" % measure, dpi = 300, bbox_inches = "tight")
    
    
    df_ld = df[df.name.isin(sig_structures_ld)]
    
    fig, ax = plt.subplots(figsize=(5,2))
    sns.stripplot(y = "name", x = measure, hue = "Condition", data = df_ld,
                orient = "h", size=7)
    # sns.boxplot(y = "name", x = "density", hue = "Condition", data = df_nm,
    #             orient = "h", showfliers=False, showcaps=False, 
    #         boxprops={'facecolor':'None'})
    plt.xlabel("Density (cells/mm$^3$)")
    plt.ylabel("Structure")
    plt.savefig("/home/wanglab/Desktop/boxplots_lower_%s.pdf" % measure, 
                dpi = 300, bbox_inches = "tight")
    
#%%
    #open ontology file
    ontology_file = "/jukebox/LightSheetData/falkner-mouse/allen_atlas/allen.json"
    with open(ontology_file) as json_file:
        ontology_dict = json.load(json_file)
        
    scale = 0.025
    flnm = "Annotated_counts"
    
    #make a list of structures you want to sum up
    struct_csv = "/jukebox/wang/zahra/cfos/structures.csv"
    dst = "/home/wanglab/Desktop/"
    
    #get list of summed structures
    sum_str = [xx[0] for xx in pd.read_csv(struct_csv,header=None).values]
    sum_str_progeny = []
    for strc in sum_str:
        progeny = []
        get_progeny(ontology_dict, strc, progeny)
        for progen in progeny:
            sum_str_progeny.append(progen)
            
    #generate data frame
    bglst=[]
    for fl in lst:
        #extract out info
        nm = os.path.basename(fl)
        #make dataframe
        df = pd.read_csv(fl+"/%s.csv" % flnm)[1:] #remove previous headers
        print(nm, df.shape)
        
        #make sure df doesnt have any summed structure names
        df = df[~df.name.isin(sum_str)]
        #sum up ontology for each structure
        for struct in df.name.values:
            progeny = []; counts = []
            get_progeny(ontology_dict, struct, progeny)
            #add structure counts to itself
            counts.append(df[df.name == struct]["counts"].values[0])
            for progen in progeny:
                try:
                    counts.append(df[df.name == progen]["counts"].values[0])
                except:
                    counts.append(0)
            df.loc[df.name == struct, "counts"] = np.sum(np.array(counts))
        #add summed structures
        for struct in sum_str:
            print(struct)
            #add summed structure to df
            progeny = []; counts = []; voxels = []
            get_progeny(ontology_dict, struct, progeny)
            for progen in progeny:
                #first counts
                try:
                    counts.append(df[df.name == progen]["counts"].values[0])
                except:
                    counts.append(0)
                try: #then voxels
                    voxels.append(df[df.name == progen]["voxels_in_structure"].values[0])
                except:
                    voxels.append(0)
            df = df.append({"name": struct, "counts": np.sum(np.array(counts)),
                            "voxels_in_structure": np.sum(np.array(voxels)), "parent_name": "root"}, 
                           ignore_index=True)
            # df.loc[df.name == struct, "counts"] = np.sum(np.array(counts))
            
        #remove summed structures progeny after you've counted them
        df = df[~df.name.isin(sum_str_progeny)]
        #remove structures you don't want to analyze
        for soi in exclude:
            df = df[df.name != soi]
            progeny = []; get_progeny(ontology_dict, soi, progeny)
            for progen in progeny:
                df = df[df.name != progen]
                
        #formatting
        df["Brain"] = nm
        df["Condition"] = conditions[nm]
        
        #append to composite dataframe
        bglst.append(df)

    df = pd.concat(bglst)
    df["counts"] = df["counts"].apply(int)

    df = df.drop(columns = ["Unnamed: 0"])
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
    
    df.to_csv(os.path.join(dst, "cell_counts_summed.csv"), index = None)
    
#%%

    df = pd.read_csv(os.path.join(dst, "cell_counts_summed.csv"), index_col = None)
    
    df_anova = pd.DataFrame()
    df_anova["name"] = df.name.unique()
    df_anova["parent_name"] = [df.loc[df.name == nm, "parent_name"].values[0]
                               for nm in df_anova["name"]]
    df_anova["tukeyhsd"] = np.ones(len(df_anova))*np.nan

    #measure to use for tests
    measure = "density"
    for nm in np.unique(df.name.values): #only gets unique names
        
        f, pval = f_oneway(df[(df.name == nm) & (df.Condition == "male-male")][measure].values, 
                     df[(df.name == nm) & (df.Condition == "female-male")][measure].values,
                     df[(df.name == nm) & (df.Condition == "male-female")][measure].values)
        
        df_anova.loc[(df_anova["name"] == nm), "anova_f"] = f
        df_anova.loc[(df_anova["name"] == nm), "anova_pval"] = pval
        
        #doing post hoc on significant structures
        if pval < 0.05:
            mc = MultiComparison(df[df.name == nm][measure].values, df[df.name == nm ].Condition.values)
            result = mc.tukeyhsd(alpha=0.05) #cutoff at 0.05
            
            df_anova.loc[(df_anova["name"] == nm), "tukeyhsd"] = np.min(result.pvalues)
                
    df_anova.to_csv(os.path.join(dst, "one_way_anova_{}.csv".format(measure)))