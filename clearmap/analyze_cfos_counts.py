#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 13:07:25 2019

@author: wanglab
"""

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
src = "/jukebox/wang/Jess/lightsheet_output/201812_development/pooled_analysis/lobuleVI_cfos/eroded_annotations_analysis"
flds = '/jukebox/wang/Jess/lightsheet_output/201812_development/forebrain/processed'
#get files
lst = [fld for fld in listdirfull(flds) if os.path.exists(os.path.join(fld, 'Annotated_counts.csv')) and fld[-4:] == "lob6" 
       or fld[-5:] == "noinj"]; lst.sort()
#conditions
nms = ['an_01_lob6',
         'an_02_lob6',
         'an_03_lob6',
         'an_04_lob6',
         'an_05_lob6',
         'an_06_lob6',
         'an_07_lob6',
         'an_09_lob6',
         'an_13_noinj',
         'an_14_noinj',
         'an_15_noinj',
         'an_16_noinj'
         ]

cond = ['Vector Control', 'DREADDs', 'Vector Control', 'DREADDs', 'Vector Control', 'DREADDs', 'DREADDs','Vector Control', 'No Injection Control',
        'No Injection Control', 'No Injection Control', 'No Injection Control']

conditions = {n:c for n,c in zip(nms, cond)}
pth = '/jukebox/wang/Jess/lightsheet_output/201812_development/pooled_analysis/lobuleVI_cfos/eroded_annotations_analysis/cell_counts_dataframe.csv'

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
        df = pd.read_csv(fl+'/Annotated_counts_eroded.csv', header = None, names = ['acronym', 
                                                                             'Count', 'Index', 'Structure', 'parent_acronym',
                                                                             'parent_name', 'Total Voxels'])[1:] #remove previous headers
        print(nm, df.shape)
        df = df.replace(np.nan, '', regex=True)
        df['Brain'] = nm
        df['Condition'] = conditions[nm]
        bglst.append(df.drop(['acronym', 'parent_acronym', 'parent_name'], axis = 1))
    
    df = pd.concat(bglst)
    df['Count'] = df['Count'].apply(int)
    df.to_csv(pth)
    
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
        structure_list = df.Structure.unique()
        brains = df.Brain.unique()
        
        print '*************'
        print c1,c2, 'len of brains: {}'.format(len(brains))
        lst = []
        for structure in structure_list:
            scount = df.loc[((df.Structure == structure) & (df.Condition == c1)), 'Count'].sum()
            smean = np.float32(df.loc[((df.Structure == structure) & (df.Condition == c1)), 'Count'].mean())
            sstd = np.float32(df.loc[((df.Structure == structure) & (df.Condition == c1)), 'Count'].std())
            ccount = df.loc[((df.Structure == structure) & (df.Condition == c2)), 'Count'].sum()
            cmean = np.float32(df.loc[((df.Structure == structure) & (df.Condition == c2)), 'Count'].mean())
            cstd = np.float32(df.loc[((df.Structure == structure) & (df.Condition == c2)), 'Count'].std())
            pval = scipy.stats.ttest_ind(np.float32(df.loc[((df.Structure == structure) & (df.Condition == c1)), 'Count']), 
                                                 np.float32(df.loc[((df.Structure == structure) & (df.Condition == c2)), 'Count']))
            try:
                mannwhit = scipy.stats.mannwhitneyu(np.float32(df.loc[((df.Structure == structure) & (df.Condition == c1)), 'Count']), 
                                                     np.float32(df.loc[((df.Structure == structure) & (df.Condition == c2)), 'Count']), alternative = 'two-sided')
            except ValueError:
                mannwhit = [0.0,1.0]
            
            wilcoxrank = scipy.stats.ranksums(np.float32(df.loc[((df.Structure == structure) & (df.Condition == c1)), 'Count']), 
                                                 np.float32(df.loc[((df.Structure == structure) & (df.Condition == c2)), 'Count']))
            lst.append((structure, scount, smean, sstd, ccount, cmean, cstd, pval[0], pval[1], mannwhit[1], wilcoxrank[1]))
            
    
        #print lst
        #make tmp df
        tdf = pd.DataFrame(data=lst, columns=['Structure', 'Stimulation count', 'Stimulation mean', 'Stimulation std', 
                                              'Control count', 'Control mean', 'Control std','tstat', 'pval', 'mannwhit', 'wilcoxrank'])
        tdf.sort_values('mannwhit')
        tdf.to_csv(os.path.join(src, 'df_with_stats_{}-{}.csv'.format(c1,c2)))
        print("saved to: {}".format(os.path.join(src, 'df_with_stats_{}-{}.csv'.format(c1,c2))))
        tdf_dct['{} vs {}'.format(c1,c2)]=tdf #this vs is important down the road
        sigs = tdf[tdf.pval<0.05].Structure.tolist()
        
    return tdf_dct, sigs

#run
tdf_dct, sigs = generate_paired_statistics(src, csv_pth)
###################################################################DONE##########################################################################################    
#%%

#helper functions
def correct_cm_to_sturctures(struc):
    '''function to correct naming issues
    '''
    if struc == 'Anterior cingulate area, ventral part, layer 6a': struc = 'Anterior cingulate area, ventral part, 6a'
    if struc == 'Anterior cingulate area, ventral part, layer 6b': struc = 'Anterior cingulate area, ventral part, 6b'
    if struc == 'Simple lobule': struc = 'Simplex lobule'
    if struc == 'Primary somatosensory area, barrel field, layer 5 ': struc = 'Primary somatosensory area, barrel field, layer 5'
    return struc

def correct_sturctures_to_cm(struc):
    '''function to correct naming issues
    '''
    if struc == 'Anterior cingulate area, ventral part, 6a': struc = 'Anterior cingulate area, ventral part, layer 6a'
    if struc == 'Anterior cingulate area, ventral part, 6b': struc = 'Anterior cingulate area, ventral part, layer 6b'
    if struc == 'Simplex lobule': struc = 'Simple lobule'
    if struc == 'Primary somatosensory area, barrel field, layer 5': struc = 'Primary somatosensory area, barrel field, layer 5 '
    return struc  

#################################################################################################################################################################


def generate_normalised_structures_list(df_pth, ann_pth, csv_pth):
    """
    generates counts normalised by volume, correct some errors in look up table
    #TODO (zahra): ask tom if this is necessary for PMA
    """
    #structures
    structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
    
    #run
    df = pd.read_csv(csv_pth)
    sois = ['Cortical subplate', 'Cortical plate', 'Cerebral nuclei', 'Thalamus', 'Hypothalamus', 'Midbrain', 'Hindbrain', 'Cerebellum', 'fiber tracts', 'ventricular systems', 'grooves']

    #add in progenitor column and add in volumes for area cell counts
    vols = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/sample_cell_count_output.xlsx')[['voxels_in_structure', 'name']]
    tdf = df.copy()
    tdf['progenitor'] = 'empty'
    tdf['Volume'] = 0.0
    scale_factor = .025 #mm/voxel
    
    for soi in sois:
        soi = [s for s in structures if s.name==soi][0]
        print soi.name
        progeny = [str(xx.name) for xx in soi.progeny]
        for progen in progeny:
            progen = correct_sturctures_to_cm(progen)
            tdf.loc[tdf['Structure']==progen,'progenitor']=soi.name
            if len(vols[vols['name']==progen]['voxels_in_structure'])>0:
                tdf.loc[tdf['Structure']==progen,'Volume']=vols[vols['name']==progen]['voxels_in_structure'].values[0]*scale_factor
    
    #drop non progen
    tdf = tdf[(tdf['progenitor'] != 'empty') & (tdf['Volume'] != 0.0)].drop(['Unnamed: 0', 'Index'], axis = 1)
    
    #add normalised column
    tdf['count_normalized_by_volume'] = tdf.apply(lambda x:x['Count']/float(x['Volume']), 1)
    
    #save both as csv and pickle for visualisation
    tdf.to_pickle(os.path.join(src, "cell_counts_dataframe_with_progenitors.p"))
    tdf.to_csv(os.path.join(src, "cell_counts_dataframe_with_progenitors.csv"))
    
    print("saved in :{}".format(src))
    
    return tdf

#run
tdf = generate_normalised_structures_list(df_pth, ann_pth, csv_pth)

###################################################################DONE##########################################################################################    
#%%

def generate_boxplots(tdf, src):
    """
    makes representative figures of cell count statistics done previously
    inputs:
        tdf: cell count dataframe (as dataframe or pickle)
        src = folder to save figs
    """
    if isinstance(tdf, basestring): tdf = pd.read_pickle(tdf)
    plt.figure(figsize=(16,8))
    g = sns.boxplot(x='progenitor', y='Count', hue='Condition', data=tdf)
    g.set_yscale('log')
    g.set_title('Cell counts boxplots by progenitor no normalization')
    plt.tight_layout()
    plt.savefig(os.path.join(src, "cell_counts_boxplots_by_progenitor_no_normalization.pdf"), dpi=300, transparent=True)
    plt.close() 

    #with cell count norm
    plt.figure(figsize=(16,8))
    g = sns.boxplot(x='progenitor', y='count_normalized_by_volume', hue='Condition', data=tdf)
    g.set_yscale('log')
    g.set_title('Cell counts boxplots by progenitor normalized by volume')
    plt.tight_layout()
    plt.savefig(os.path.join(src, "cell_counts_boxplots_by_progenitor_normalized_by_volume.pdf"), dpi=300, transparent=True)
    plt.close()

    #total counts per brain
    plt.figure(figsize=(16,8))
    g=sns.pairplot(tdf, hue='Brain')
    g.set(xscale='log', yscale="log")
    #plt.tight_layout()
    plt.savefig(os.path.join(src, "pairplot_by_brain.pdf"), dpi=300, transparent=True)
    plt.close()

    plt.figure(figsize=(16,8))
    g=sns.pairplot(tdf, hue='Condition')
    g.set(xscale='log', yscale="log")
    #plt.tight_layout()
    plt.savefig(os.path.join(src, "pairplot_by_condition.pdf"), dpi=300, transparent=True)
    plt.close()

    #sum
    #g=tdf.groupby('Condition').sum().plot(kind='bar')
    #g.set_yscale('log')
    plt.figure(figsize=(16,8))
    g = sns.boxplot(x='Brain', y='Count', data=tdf)
    g.set_yscale('log')
    g.set_xticklabels(g.axes.get_xticklabels(), rotation=30)
    g.set_title('Structure total counts per brain')
    plt.tight_layout()
    plt.savefig(os.path.join(src, "sum_by_brain.pdf"), dpi=300, transparent=True)
    plt.close()
    
    print("saved in : {}".format(src))
    
#run
generate_boxplots(tdf, src)

###################################################################DONE##########################################################################################        
#%%
#since you can't really to stats on n=2. Possibly make ratios of structure in condition a with b. Then pool. 
#Look for bi/tri-modality as this virus might have brain areas that it better infects than others.
def normalise_counts_by(src, tdf):
    """ normalise counts by each condition to see bi/tri-modality """
    
    tdf['condition_structure_mean'] = 0.0
    ndf = pd.DataFrame(index = tdf['Structure'].unique(), columns = tdf['Condition'].unique()) #values = condition_structure_mean
    ndf[:] = np.nan
    for s in tdf.Structure.unique():
        for cond in tdf['Condition'].unique():
            mn = tdf[(tdf['Condition']==cond) &(tdf['Structure']==s)]['Count'].mean()
            tdf.loc[(tdf['Condition']==cond) &(tdf['Structure']==s), 'structure_mean'] = mn
            ndf.loc[ndf.index==s, ndf.columns==cond] = mn
    
    
    for cond in tdf['Condition'].unique():
        nndf=ndf.copy()
        others = tdf['Condition'].unique()
        others = np.delete(others, np.argwhere(others==cond))
        for o in others:
            print o, cond
            nndf['{}_{}'.format(o,cond)] = nndf.apply(lambda x: x[o]/float(x[cond]), 1)
        nndf = nndf[nndf.columns.difference(ndf.columns)]
        #deal with Nan's in structs no represented
        fig = plt.figure(figsize=(16,8))
        g=sns.boxplot(data=nndf)    
        #g=sns.swarmplot(data=nndf)    
        g.axes.axhline(y=1, xmin=0.0, xmax=1.0, color='k', linestyle='--')
        g.set_yscale('log')
        plt.title('Counts normalized by {}.pdf'.format(cond))
        plt.tight_layout()
        plt.savefig(os.path.join(src, "boxplot_counts_normalized_by_{}.pdf".format(cond)), dpi=300, transparent=True)
        plt.close()
        
        #density
        fig = plt.figure(figsize=(16,8))
        fig.add_subplot(2,1,1)
        for c in nndf.columns:
            g=sns.kdeplot(nndf[c].dropna(), shade=True)
        g.axes.set_xlabel('Counts')
        g.axes.set_ylabel('Ratio of counts')
        plt.title('Counts normalized by {}.pdf'.format(cond))
        g.axes.set_xlim(0,g.axes.get_xlim()[1])
        plt.tight_layout()
        fig.add_subplot(2,1,2)
        for c in nndf.columns:
            g=sns.kdeplot(nndf[c].dropna(), shade=True)
        g.axes.set_ylim(g.axes.get_ylim()[0]/4, g.axes.get_ylim()[1]/4)
        g.axes.set_xlim(g.axes.get_xlim()[0]/4, g.axes.get_xlim()[1]/4)
        g.axes.set_xlabel('Counts')
        g.axes.set_ylabel('Ratio of counts')
        plt.title('Counts normalized by {} zoomed.pdf'.format(cond))
        plt.tight_layout()
        plt.savefig(os.path.join(src, "kde_counts_normalized_by_{}.pdf".format(cond)), dpi=300, transparent=True)
        plt.close()
        
        print("saved in : {}".format(src))

#run
normalise_counts_by(src, tdf)    

###################################################################DONE##########################################################################################        
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
out = '/jukebox/wang/Jess/lightsheet_output/201812_development/pooled_analysis/images'
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
    groupC = [fld for fld in listdirfull(flds) if fld[-4:] == "lob6" and conditions[os.path.basename(fld)] == "DREADDs"]; groupC.sort()
    groupD = [fld for fld in listdirfull(flds) if fld[-4:] == "lob6" and conditions[os.path.basename(fld)] == "Vector Control"]; groupD.sort()
    groupE = [fld for fld in listdirfull(flds) if fld[-5:] == "noinj"]; groupE.sort() 
    
    group_c = [xx+'/cells_heatmap.tif' for xx in groupC]  
    group_d = [xx+'/cells_heatmap.tif' for xx in groupD]  
    group_e = [xx+'/cells_heatmap.tif' for xx in groupE]  
    
    grp_c = stat.readDataGroup(group_c)
    grp_d = stat.readDataGroup(group_d)
    grp_e = stat.readDataGroup(group_e)
    
    #Generated average and standard deviation maps
    ##############################################
    grp_ca = np.mean(grp_c, axis = 0)
    grp_cs = np.std(grp_c, axis = 0)
    
    grp_da = np.mean(grp_d, axis = 0)
    grp_ds = np.std(grp_d, axis = 0)
    
    grp_ea = np.mean(grp_e, axis = 0)
    grp_es = np.std(grp_e, axis = 0)
    
    io.writeData(os.path.join(src, 'group_c_mean.raw'), rsp.sagittalToCoronalData(grp_ca))
    io.writeData(os.path.join(src, 'group_c_std.raw'), rsp.sagittalToCoronalData(grp_cs))
    
    io.writeData(os.path.join(src, 'group_d_mean.raw'), rsp.sagittalToCoronalData(grp_da))
    io.writeData(os.path.join(src, 'group_d_std.raw'), rsp.sagittalToCoronalData(grp_ds))
    
    io.writeData(os.path.join(src, 'group_e_mean.raw'), rsp.sagittalToCoronalData(grp_ea))
    io.writeData(os.path.join(src, 'group_e_std.raw'), rsp.sagittalToCoronalData(grp_es))
    
    #Generate the p-values map
    ##########################
    #first comparison
    #pcutoff: only display pixels below this level of significance
    pvals, psign = stat.tTestVoxelization(grp_c.astype('float'), grp_d.astype('float'), signed = True, pcutoff = 0.05)
    
    #color the p-values according to their sign (defined by the sign of the difference of the means between the 2 groups)
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    io.writeData(os.path.join(src, 'pvalues_groupCvsD.tif'), rsp.sagittalToCoronalData(pvalsc.astype('float32')));
    
    #second comparison
    pvals, psign = stat.tTestVoxelization(grp_c.astype('float'), grp_e.astype('float'), signed = True, pcutoff = 0.05)
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    io.writeData(os.path.join(src, 'pvalues_groupCvsE.tif'), rsp.sagittalToCoronalData(pvalsc.astype('float32')))

#run
generate_p_value_maps(src)
###################################################################DONE##########################################################################################        

