#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 18:28:08 2019

@author: wanglab
"""

src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/transformed_points"
lob6 = [ '20180608_jg75', '20180410_jg48_bl6_lob6a_01', '20170115_tp_bl6_lob6a_rpv_03', '20170115_tp_bl6_lob6b_ml_04',
       '20180410_jg50_bl6_lob6b_03', '20170115_tp_bl6_lob6a_1000r_02','20170115_tp_bl6_lob6a_500r_01' ]
post_transformed = [os.path.join(src, os.path.join(xx, "transformed_points/posttransformed_zyx_voxels.npy")) for xx in lob6]

def transformed_cells_to_allen(fld, ann, dst, fl_nm):
    """ consolidating to one function bc then no need to copy/paste """
    dct = {}
    
    for fl in fld:
        brain = os.path.basename(os.path.dirname(os.path.dirname(fl)))
        point_lst = transformed_pnts_to_allen_helper_func(np.load(fl), tifffile.imread(ann), order = "ZYX")
        df = count_structure_lister(id_table, *point_lst).fillna(0)
        #for some reason duplicating columns, so use this
        nm_cnt = pd.Series(df.cell_count.values, df.name.values).to_dict()
        fl_name = brain
        dct[fl_name]= nm_cnt
        
    #unpack
    index = dct[list(dct.keys())[0]].keys()
    columns = dct.keys()
    data = np.asarray([[dct[col][idx] for idx in index] for col in columns])
    df = pd.DataFrame(data.T, columns=columns, index=index)
    
    #save before adding projeny counts at each level
    df.to_pickle(os.path.join(dst, fl_nm))
    
    return os.path.join(dst, fl_nm)

ann_pth = "/home/wanglab/mounts/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_60um_erosion.tif"
dst = "/home/wanglab/Desktop"
pth = transformed_cells_to_allen(post_transformed, ann_pth, dst, "tracing_nc_counts.p")

#%%
def get_cell_n_density_counts(brains, structure, structures, cells_regions, scale_factor = 0.025):
    """ consolidating to one function bc then no need to copy/paste """
    #get cell counts adn densities
    #get densities for all the structures
    df = pd.read_excel("/home/wanglab/mounts/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", index_col = None)
    df = df.drop(columns = ["Unnamed: 0"])
    df = df.sort_values(by = ["name"])
    
    #make new dict - for all brains
    cells_pooled_regions = {} #for raw counts
    volume_pooled_regions = {} #for density
    
    for brain in brains:    
        #make new dict - this is for EACH BRAIN
        c_pooled_regions = {}
        d_pooled_regions = {}
        
        for soi in structure:
            try:
                soi = [s for s in structures if s.name==soi][0]
                counts = [] #store counts in this list
                volume = [] #store volume in this list
                for k, v in cells_regions[brain].items():
                    if k == soi.name:
                        counts.append(v)
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))
                progeny = [str(xx.name) for xx in soi.progeny]
                #now sum up progeny
                if len(progeny) > 0:
                    for progen in progeny:
                        for k, v in cells_regions[brain].items():
                            if k == progen and progen != "Primary somatosensory area, unassigned, layer 4,5,6":
                                counts.append(v)
                                #add to volume list from LUT
                                volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0])
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
            except:
                for k, v in cells_regions[brain].items():
                    if k == soi:
                        counts.append(v)                    
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))                c_pooled_regions[soi] = np.sum(np.asarray(counts))
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
                        
        #add to big dict
        cells_pooled_regions[brain] = c_pooled_regions
        volume_pooled_regions[brain] = d_pooled_regions
    #making the proper array per brain where regions are removed
    cell_counts_per_brain = []
    #initialise dummy var
    i = []
    for k,v in cells_pooled_regions.items():
        dct = cells_pooled_regions[k]
        for j,l in dct.items():
            i.append(l)  
        cell_counts_per_brain.append(i)
        #re-initialise for next
        i = []  
    cell_counts_per_brain = np.asarray(cell_counts_per_brain)
    
    volume_per_brain = []
    #initialise dummy var
    i = []
    for k,v in volume_pooled_regions.items():
        dct = volume_pooled_regions[k]
        for j,l in dct.items():
            i.append(l)  
        volume_per_brain.append(i)
        #re-initialise for next
        i = []  
    volume_per_brain = np.asarray(volume_per_brain)*(scale_factor**3)
    #calculate denisty
    density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])
    
    return cell_counts_per_brain, density_per_brain

#making dictionary of cells by region
cells_regions = pckl.load(open("/home/wanglab/Desktop/tracing_nc_counts.p", "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

nc_areas = ["Central amygdalar nucleus", "Intercalated amygdalar nucleus", "Anterior amygdalar area",
            "Medial amygdalar nucleus", "Basolateral amygdalar nucleus", "Lateral amygdalar nucleus",
            "Basomedial amygdalar nucleus", "Posterior amygdalar nucleus"]

cell_counts_per_brain, density_per_brain = get_cell_n_density_counts(lob6, nc_areas, structures, cells_regions, 
                                                                     scale_factor = 0.020)

#%%
df = pd.DataFrame(cell_counts_per_brain)
df.columns = nc_areas
df.index = lob6
df.to_csv("/home/wanglab/Desktop/lob6_tracing_cell_counts_amyg.csv")

df = pd.DataFrame(density_per_brain)
df.columns = nc_areas
df.index = lob6
df.to_csv("/home/wanglab/Desktop/lob6_tracing_density_amyg.csv")
