#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 18:17:23 2019

@author: wanglab
"""

from __future__ import division
import os, numpy as np
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from skimage.external import tifffile
import seaborn as sns, pandas as pd, matplotlib.pyplot as plt
import scipy, itertools
from skimage.exposure import equalize_hist, adjust_gamma
from tools.utils.io import listdirfull
from tools.analysis.network_analysis import make_structure_objects
from tools.utils.overlay import tile
os.chdir("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.cluster.utils import load_kwargs
import ClearMap.IO.IO as io
import ClearMap.Analysis.Statistics as stat
import ClearMap.Alignment.Resampling as rsp

sns.set_style('white')

#make inputs
src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis"
flds = '/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/processed'
#get files
lst = [fld for fld in listdirfull(flds) if os.path.exists(os.path.join(fld, 'Annotated_counts.csv'))]; lst.sort()
#conditions
nms = ['an1',
         'an2',
         'an3',
         'an4',
         'an5',
         'an6',
         'an7',
         'an8',
         'an9',
         'an10',
         'an11',
         'an12',
         'an13',
         'an14',
         'an15',
         'an16',
         'an18',
         'an19',
         'an20',
         'an21',
         'an22',
         'an23',
         'an24',
         'an25',
         'an26',
         'an27',
         'an28',
         'an29',
         'an30',
         'an31',
         'an32',
         'an33'
     ]

cond = ['CNO_control_no_reversal', 'CNO_control_no_reversal', 'CNO_control_no_reversal', 
        'CNO_control_no_reversal', 'CNO_control_no_reversal', 'CNO_control_no_reversal',
        'CNO_control_no_reversal', 'CNO_control_no_reversal', 
        "CNO_control_reversal", "CNO_control_reversal", "CNO_control_reversal", "CNO_control_reversal",
        "CNO_control_reversal", "CNO_control_reversal", "CNO_control_reversal", "CNO_control_reversal",
        "CNO_control_reversal", "DREADDs",
        "DREADDs", "DREADDs", "DREADDs", "DREADDs", "DREADDs", "DREADDs", "DREADDs", "DREADDs", "DREADDs", 
        "homecage_control", "homecage_control", "homecage_control", "homecage_control", "homecage_control"]

conditions = {n:c for n,c in zip(nms, cond)}
pth = os.path.join(src, 'cell_counts_dataframe.csv')

df_pth = '/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx'
ann_pth = '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_60um_erosion.tif'
atl_pth = '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif'

#%%

def generate_data_frame(conditions, lst, pth):
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
        df = pd.read_csv(fl+'/Annotated_counts.csv')[1:] #remove previous headers
        print(nm, df.shape)
        df = df.replace(np.nan, '', regex=True)
        df['Brain'] = nm
        df['Condition'] = conditions[nm]
        bglst.append(df.drop(['acronym', 'parent_acronym', 'parent_name'], axis = 1))
    
    df = pd.concat(bglst)
    df['counts'] = df['counts'].apply(int)
    df.drop(columns = ["Unnamed: 0"]).to_csv(pth, index = None)
    
    return pth

#run
csv_pth = generate_data_frame(conditions, lst, pth)
###################################################################DONE##########################################################################################
#%%

def generate_paired_statistics(src, csv_pth):
    """
    generates paried t tests, and Mann Whitney and Wilcox Rank test results from pooled csv counts
    returns:
        tdf_dct: data frame with comparison tests
        sigs: significant structures (p < 0.05)
    """
    
    df = pd.read_csv(csv_pth)
    tdf_dct={}
    
    for c1,c2 in itertools.combinations(df['Condition'].unique(), 2):    
        df = pd.read_csv(csv_pth)
        df = df[df['Condition'].isin([c1,c2])]                                                                                                             
        structure_list = df.name.unique()
        brains = df.Brain.unique()
        
        print '*************'
        print c1,c2, 'len of brains: {}'.format(len(brains))
        lst = []
        for structure in structure_list:
            scount = df.loc[((df.name == structure) & (df.Condition == c1)), 'counts'].sum()
            smean = np.float32(df.loc[((df.name == structure) & (df.Condition == c1)), 'counts'].mean())
            sstd = np.float32(df.loc[((df.name == structure) & (df.Condition == c1)), 'counts'].std())
            ccount = df.loc[((df.name == structure) & (df.Condition == c2)), 'counts'].sum()
            cmean = np.float32(df.loc[((df.name == structure) & (df.Condition == c2)), 'counts'].mean())
            cstd = np.float32(df.loc[((df.name == structure) & (df.Condition == c2)), 'counts'].std())
            pval = scipy.stats.ttest_ind(np.float32(df.loc[((df.name == structure) & (df.Condition == c1)), 'counts']), 
                                                 np.float32(df.loc[((df.name == structure) & (df.Condition == c2)), 'counts']))
            try:
                mannwhit = scipy.stats.mannwhitneyu(np.float32(df.loc[((df.name == structure) & (df.Condition == c1)), 'counts']), 
                                                     np.float32(df.loc[((df.name == structure) & (df.Condition == c2)), 'counts']), alternative = 'two-sided')
            except ValueError:
                mannwhit = [0.0,1.0]
            
            wilcoxrank = scipy.stats.ranksums(np.float32(df.loc[((df.name == structure) & (df.Condition == c1)), 'counts']), 
                                                 np.float32(df.loc[((df.name == structure) & (df.Condition == c2)), 'counts']))
            lst.append((structure, scount, smean, sstd, ccount, cmean, cstd, pval[0], pval[1], mannwhit[1], wilcoxrank[1]))
            
    
        #print lst
        #make tmp df
        tdf = pd.DataFrame(data=lst, columns=['name', 'Stimulation count', 'Stimulation mean', 'Stimulation std', 
                                              'Control count', 'Control mean', 'Control std','tstat', 'pval', 'mannwhit', 'wilcoxrank'])
        tdf.sort_values('mannwhit')
        tdf.to_csv(os.path.join(src, 'df_with_stats_{}-{}.csv'.format(c1,c2)))
        print("saved to: {}".format(os.path.join(src, 'df_with_stats_{}-{}.csv'.format(c1,c2))))
        tdf_dct['{} vs {}'.format(c1,c2)]=tdf #this vs is important down the road
        sigs = tdf[tdf.pval<0.05].name.tolist()
        
    return tdf_dct, sigs

#run
tdf_dct, sigs = generate_paired_statistics(src, csv_pth)
###################################################################DONE##########################################################################################    
#%%
#
##helper functions
#def correct_cm_to_sturctures(struc):
#    '''function to correct naming issues
#    '''
#    if struc == 'Anterior cingulate area, ventral part, layer 6a': struc = 'Anterior cingulate area, ventral part, 6a'
#    if struc == 'Anterior cingulate area, ventral part, layer 6b': struc = 'Anterior cingulate area, ventral part, 6b'
#    if struc == 'Simple lobule': struc = 'Simplex lobule'
#    if struc == 'Primary somatosensory area, barrel field, layer 5 ': struc = 'Primary somatosensory area, barrel field, layer 5'
#    return struc
#
#def correct_sturctures_to_cm(struc):
#    '''function to correct naming issues
#    '''
#    if struc == 'Anterior cingulate area, ventral part, 6a': struc = 'Anterior cingulate area, ventral part, layer 6a'
#    if struc == 'Anterior cingulate area, ventral part, 6b': struc = 'Anterior cingulate area, ventral part, layer 6b'
#    if struc == 'Simplex lobule': struc = 'Simple lobule'
#    if struc == 'Primary somatosensory area, barrel field, layer 5': struc = 'Primary somatosensory area, barrel field, layer 5 '
#    return struc  
#
#
#
#def generate_normalised_structures_list(df_pth, ann_pth, csv_pth):
#    """
#    generates counts normalised by volume, correct some errors in look up table
#    #TODO (zahra): ask tom if this is necessary for PMA
#    """
#    #structures
#    structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
#    
#    #run
#    df = pd.read_csv(csv_pth)
#    sois = ['Cortical subplate', 'Cortical plate', 'Cerebral nuclei', 'Thalamus', 'Hypothalamus', 'Midbrain', 'Hindbrain', 'Cerebellum', 'fiber tracts', 'ventricular systems', 'grooves']
#
#    #add in progenitor column and add in volumes for area cell counts
#    vols = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/sample_cell_count_output.xlsx')[['voxels_in_structure', 'name']]
#    tdf = df.copy()
#    tdf['progenitor'] = 'empty'
#    tdf['Volume'] = 0.0
#    scale_factor = .025 #mm/voxel
#    
#    for soi in sois:
#        soi = [s for s in structures if s.name==soi][0]
#        print soi.name
#        progeny = [str(xx.name) for xx in soi.progeny]
#        for progen in progeny:
#            progen = correct_sturctures_to_cm(progen)
#            tdf.loc[tdf['name']==progen,'progenitor']=soi.name
#            if len(vols[vols['name']==progen]['voxels_in_structure'])>0:
#                tdf.loc[tdf['name']==progen,'Volume']=vols[vols['name']==progen]['voxels_in_structure'].values[0]*scale_factor
#    
#    #drop non progen
#    tdf = tdf[(tdf['progenitor'] != 'empty') & (tdf['Volume'] != 0.0)].drop(['Unnamed: 0', 'Index'], axis = 1)
#    
#    #add normalised column
#    tdf['count_normalized_by_volume'] = tdf.apply(lambda x:x['Count']/float(x['Volume']), 1)
#    
#    #save both as csv and pickle for visualisation
#    tdf.to_pickle(os.path.join(src, "cell_counts_dataframe_with_progenitors.p"))
#    tdf.to_csv(os.path.join(src, "cell_counts_dataframe_with_progenitors.csv"))
#    
#    print("saved in :{}".format(src))
#    
#    return tdf
#
##run
#tdf = generate_normalised_structures_list(df_pth, ann_pth, csv_pth)
#
##%%
#
#def generate_boxplots(tdf, src):
#    """
#    makes representative figures of cell count statistics done previously
#    inputs:
#        tdf: cell count dataframe (as dataframe or pickle)
#        src = folder to save figs
#    """
#    if isinstance(tdf, basestring): tdf = pd.read_pickle(tdf)
#    plt.figure(figsize=(16,8))
#    g = sns.boxplot(x='progenitor', y='Count', hue='Condition', data=tdf)
#    g.set_yscale('log')
#    g.set_title('Cell counts boxplots by progenitor no normalization')
#    plt.tight_layout()
#    plt.savefig(os.path.join(src, "cell_counts_boxplots_by_progenitor_no_normalization.pdf"), dpi=300, transparent=True)
#    plt.close() 
#
#    #with cell count norm
#    plt.figure(figsize=(16,8))
#    g = sns.boxplot(x='progenitor', y='count_normalized_by_volume', hue='Condition', data=tdf)
#    g.set_yscale('log')
#    g.set_title('Cell counts boxplots by progenitor normalized by volume')
#    plt.tight_layout()
#    plt.savefig(os.path.join(src, "cell_counts_boxplots_by_progenitor_normalized_by_volume.pdf"), dpi=300, transparent=True)
#    plt.close()
#
#    #total counts per brain
#    plt.figure(figsize=(16,8))
#    g=sns.pairplot(tdf, hue='Brain')
#    g.set(xscale='log', yscale="log")
#    #plt.tight_layout()
#    plt.savefig(os.path.join(src, "pairplot_by_brain.pdf"), dpi=300, transparent=True)
#    plt.close()
#
#    plt.figure(figsize=(16,8))
#    g=sns.pairplot(tdf, hue='Condition')
#    g.set(xscale='log', yscale="log")
#    #plt.tight_layout()
#    plt.savefig(os.path.join(src, "pairplot_by_condition.pdf"), dpi=300, transparent=True)
#    plt.close()
#
#    #sum
#    #g=tdf.groupby('Condition').sum().plot(kind='bar')
#    #g.set_yscale('log')
#    plt.figure(figsize=(16,8))
#    g = sns.boxplot(x='Brain', y='Count', data=tdf)
#    g.set_yscale('log')
#    g.set_xticklabels(g.axes.get_xticklabels(), rotation=30)
#    g.set_title('Structure total counts per brain')
#    plt.tight_layout()
#    plt.savefig(os.path.join(src, "sum_by_brain.pdf"), dpi=300, transparent=True)
#    plt.close()
#    
#    print("saved in : {}".format(src))
#    
##run
#generate_boxplots(tdf, src)

#%% #look at cross sections

def check_registration_cross_sections(out):
    for z in [100,200,300,400,500]:
        print z
        nm_im = {}
        for fld in lst:
            kwargs = load_kwargs(fld)
            vol = [xx for xx in kwargs['volumes'] if xx.ch_type == 'cellch'][0]
            fl = [fl for fl in listdirfull(vol.full_sizedatafld_vol) if str(z).zfill(4) in fl][0]
            nm_im[os.path.basename(fld)] = fl
            
        dst = os.path.join(out, 'cell_ch_z{}.png'.format(str(z).zfill(4)))
        tile(src = [adjust_gamma(tifffile.imread(xx), gamma=.6,gain=3) for xx in nm_im.values()], subtitles=[xx for xx in nm_im.keys()], dst = dst)
        
    
    #check reg
    nm_im = {}
    for fld in lst:
        kwargs = load_kwargs(fld)
        fl = os.path.join(fld, 'clearmap_cluster_output', 'elastix_auto_to_atlas', 'result.1.tif')
        if os.path.exists(fl):
            nm_im[os.path.basename(fld)] = fl
            
    #read once
    ims = [equalize_hist(tifffile.imread(xx)) for xx in nm_im.values()]
    for z in [50,100,150,200,250,300,350,400]:
        print z
        dst = os.path.join(out, 'regqc_z{}.png'.format(str(z).zfill(4)))
        tile(src = [i[z] for i in ims], subtitles=[xx for xx in nm_im.keys()], dst = dst)
        
    print("saved in : {}".format(out))
        
#run
out = os.path.join(src, "images")
if not os.path.exists(out): os.mkdir(out)
check_registration_cross_sections(out)
###################################################################DONE##########################################################################################        
#%%
#pooled analysis to make p-value maps 

def generate_p_value_maps(src):
    
    """ 
    generates p-value maps as per ClearMap/analysis.py
    #TODO: generalise function
    """
    #Load the data (heat maps generated previously )
    #make groups
    groupA = [fld for fld in listdirfull(flds) if conditions[os.path.basename(fld)] == "homecage_control"]; groupA.sort() 
    groupB = [fld for fld in listdirfull(flds) if conditions[os.path.basename(fld)] == "CNO_control_no_reversal"]; groupB.sort()
    groupC = [fld for fld in listdirfull(flds) if conditions[os.path.basename(fld)] == "CNO_control_reversal"]; groupC.sort()
    groupD = [fld for fld in listdirfull(flds) if conditions[os.path.basename(fld)] == "DREADDs"]; groupC.sort()

    group_a = [xx+'/cells_heatmap.tif' for xx in groupA]
    group_b = [xx+'/cells_heatmap.tif' for xx in groupB]  
    group_c = [xx+'/cells_heatmap.tif' for xx in groupC]  
    group_d = [xx+'/cells_heatmap.tif' for xx in groupD]  
    
    
    grp_a = stat.readDataGroup(group_a)
    grp_b = stat.readDataGroup(group_b)
    grp_c = stat.readDataGroup(group_c)
    grp_d = stat.readDataGroup(group_d)
    
    #Generated average and standard deviation maps
    ##############################################
    grp_aa = np.mean(grp_a, axis = 0)
    grp_as = np.std(grp_a, axis = 0)
    
    grp_ba = np.mean(grp_b, axis = 0)
    grp_bs = np.std(grp_b, axis = 0)
    
    grp_ca = np.mean(grp_c, axis = 0)
    grp_cs = np.std(grp_c, axis = 0)
    
    grp_da = np.mean(grp_d, axis = 0)
    grp_ds = np.std(grp_d, axis = 0)
    
    io.writeData(os.path.join(src, 'group_a_mean.raw'), rsp.sagittalToCoronalData(grp_aa))
    io.writeData(os.path.join(src, 'group_a_std.raw'), rsp.sagittalToCoronalData(grp_as))
    
    io.writeData(os.path.join(src, 'group_b_mean.raw'), rsp.sagittalToCoronalData(grp_ba))
    io.writeData(os.path.join(src, 'group_b_std.raw'), rsp.sagittalToCoronalData(grp_bs))
    
    io.writeData(os.path.join(src, 'group_c_mean.raw'), rsp.sagittalToCoronalData(grp_ca))
    io.writeData(os.path.join(src, 'group_c_std.raw'), rsp.sagittalToCoronalData(grp_cs))
    
    io.writeData(os.path.join(src, 'group_d_mean.raw'), rsp.sagittalToCoronalData(grp_da))
    io.writeData(os.path.join(src, 'group_d_std.raw'), rsp.sagittalToCoronalData(grp_ds))
    
    #Generate the p-values map
    ##########################
    #first comparison
    #pcutoff: only display pixels below this level of significance
    pvals, psign = stat.tTestVoxelization(grp_a.astype('float'), grp_d.astype('float'), signed = True, pcutoff = 0.05)
    
    #color the p-values according to their sign (defined by the sign of the difference of the means between the 2 groups)
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    io.writeData(os.path.join(src, 'pvalues_homecage_control_vs_DREADDs.tif'), rsp.sagittalToCoronalData(pvalsc.astype('float32')));
    
    #second comparison
    pvals, psign = stat.tTestVoxelization(grp_a.astype('float'), grp_b.astype('float'), signed = True, pcutoff = 0.05)
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    io.writeData(os.path.join(src, 'pvalues_homecage_control_vs_CNO_control_no_reversal.tif'), rsp.sagittalToCoronalData(pvalsc.astype('float32')))
    
    #third comparison
    pvals, psign = stat.tTestVoxelization(grp_b.astype('float'), grp_c.astype('float'), signed = True, pcutoff = 0.05)
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    io.writeData(os.path.join(src, 'pvalues_CNO_control_no_reversal_vs_CNO_control_reversal.tif'), rsp.sagittalToCoronalData(pvalsc.astype('float32')))
    
    #fourth comparison
    pvals, psign = stat.tTestVoxelization(grp_c.astype('float'), grp_d.astype('float'), signed = True, pcutoff = 0.05)
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    io.writeData(os.path.join(src, 'pvalues_CNO_control_reversal_vs_DREADDs.tif'), rsp.sagittalToCoronalData(pvalsc.astype('float32')))

#run
src = os.path.join(src, "p_value_maps")
if not os.path.exists(src): os.mkdir(src)
generate_p_value_maps(src)

###################################################################DONE##########################################################################################        
#%%

def generate_percent_counts_and_density_per_region(src, csv_pth):
    
    """ generates another column in the dataframe that is just # counts / total counts in brain """
    
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
    
    df["percent_per_total_counts (%)"] = df_percent
    df["volume"] = df[df.voxels_in_structure > 0].apply(lambda x: x.voxels_in_structure*(scale**3), 1)
    df["density (cells/mm3)"] = df[df.voxels_in_structure > 0].apply(lambda x:x.counts/(float(x.voxels_in_structure*(scale**3))), 1)
    #export
    df.to_csv(os.path.join(src, "cell_counts_dataframe_w_percents_density.csv"))
    
    return os.path.join(src, "cell_counts_dataframe_w_percents_density.csv")

#run
percent_density_csv_pth = generate_percent_counts_and_density_per_region(src, csv_pth)
    
###################################################################DONE##########################################################################################        
#%%

def pool_regions(sois, df_pth, ann_pth, src, csv_pth):

    """ pools regions together based on progenitors """    
    
    #build structures class
    structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
    
    #set variables
    df = pd.read_csv(csv_pth)
    

    #add in progenitor column
    tdf = df.copy()
    tdf["parent"] = "empty"
        
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print(soi.name)
        progeny = [str(xx.name) for xx in soi.progeny]
        if len(progeny) > 0:
            for progen in progeny:
                tdf.loc[tdf.name == progen,"parent"]=soi.name
    
    #drop non progen
    tdf = tdf[(tdf["parent"] != "empty") & (tdf["volume"] != 0.0)]
    
 
    #drop unnamed
    tdf = tdf.drop(columns = ["Unnamed: 0"])
    #save both as csv and pickle for visualisation
    tdf.to_pickle(os.path.join(src, "cell_counts_dataframe_with_progenitors.p"))
    tdf.to_csv(os.path.join(src, "cell_counts_dataframe_with_progenitors.csv"), header = True, index = None)
    
    #aggregate counts
    df_new = tdf.groupby(["parent", "Brain", "Condition"])['counts'].sum()
    df_new = pd.DataFrame(df_new)
    df_new['percent_per_total_counts (%)'] = tdf.groupby(["parent", "Brain"])['percent_per_total_counts (%)'].sum()
    df_new['density (cells/mm3)'] = tdf.groupby(["parent", "Brain"])['density (cells/mm3)'].sum()
    
    #save as new csv 
    df_new.to_csv(os.path.join(src, "summed_parents_cell_counts_dataframe.csv"), header = True)


    print("saved in :{}".format(src))
    
    return tdf

#give list of structures you want to pool
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area", "Entorhinal area", "Thalamus, sensory-motor cortex related", 
            "Thalamus, polymodal association cortex related", "Periventricular zone", "Periventricular region", "Hypothalamic medial zone", 
            "Hypothalamic lateral zone", "Midbrain, sensory related", "Midbrain, motor related", "Midbrain, behavioral state related", 
            "Hippocampal formation", "Striatum", "Pallidum", "Periaqueductal gray", "Subiculum", "Main olfactory bulb", "Olfactory areas"
            ]

#run
pool_regions(sois, df_pth, ann_pth, src, percent_density_csv_pth)