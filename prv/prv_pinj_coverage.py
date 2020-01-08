#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 13:32:31 2020

@author: wanglab
"""
#get percent coverage of each
for injdst in ['/home/wanglab/wang/pisano/tracing_output/antero_4x_analysis/201902_injection/nc', '/home/wanglab/wang/pisano/tracing_output/antero_4x_analysis/201902_injection/thalamus']:
    print injdst
    #get percentages
    #cb = [s for s in structures if 'Cerebellar nuclei' in s.name][0]
    #cb = [s for s in structures if 'Cerebellar cortex' in s.name][0]
    df = pd.read_csv(os.path.join(injdst, 'voxel_counts.csv'))
    cb = [s for s in structures if 'Cerebellum' in s.name][0]
    verm_order = [u'Lingula (I)',u'Lobule II',u'Lobule III',u'Lobules IV-V',u'Declive (VI), subdivision A',u'Declive (VI), subdivision B',u'Folium-tuber vermis (VII)',u'Pyramus (VIII)'u'Uvula (IX)',u'Nodulus (X)']
    hemi_order = [u'Simplex lobule',u'Crus 1',u'Crus 2',u'Paramedian lobule',u'Copula pyramidis']#,u'Paraflocculus',u'Flocculus']
    dcn_order = [u'Fastigial nucleus',u'Interposed nucleus',u'Dentate nucleus']#,u'Vestibulocerebellar nucleus']
    #crop...orientation is done AFTER saving individual ones...
    ann = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif')[:,423:,:]
    atl = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif')[:,423:,:]
    number_of_voxels = 1
    all_arrs = np.asarray([tifffile.imread(fl).astype('bool').astype('uint8') for fl in listall(injdst+'_arrays', '.tif')])
    sm = all_arrs.sum(0)
    assert sm.shape == ann.shape
    cbs = {xx.name: xx for xx in cb.progeny if xx.voxels_in_structure>0}
    data=[]
    
    for xx in cbs.values():
        print xx.name
        s_mask = np.where(ann==xx.id)
        vox_in_structure = s_mask[0].shape[0]
        _0 = np.where(sm[s_mask]==0)[0].shape[0]
        _1 = np.where(sm[s_mask]>=1)[0].shape[0]
        _2 = np.where(sm[s_mask]>=2)[0].shape[0]
        _3 = np.where(sm[s_mask]>=3)[0].shape[0]
        _4 = np.where(sm[s_mask]>=4)[0].shape[0]
        data.append([xx.name, 'both', vox_in_structure, _0/float(vox_in_structure), _1/float(vox_in_structure), _2/float(vox_in_structure), _3/float(vox_in_structure), _4/float(vox_in_structure)])
        if xx.name == 'Crus 1': print _1
        #only for hemis
        if xx.name in hemi_order+dcn_order:
            #right
            rann = ann[:270]
            rsm = sm[:270]
            s_mask = np.where(rann==xx.id)
            vox_in_structure = s_mask[0].shape[0]
            _0 = np.where(rsm[s_mask]==0)[0].shape[0]
            _1 = np.where(rsm[s_mask]>=1)[0].shape[0]
            _2 = np.where(rsm[s_mask]>=2)[0].shape[0]
            _3 = np.where(rsm[s_mask]>=3)[0].shape[0]
            _4 = np.where(rsm[s_mask]>=4)[0].shape[0]
            data.append([xx.name, 'right',vox_in_structure, _0/float(vox_in_structure), _1/float(vox_in_structure), _2/float(vox_in_structure), _3/float(vox_in_structure), _4/float(vox_in_structure)])
            if xx.name == 'Crus 1': print _1
            
            #left
            lann = ann[270:]
            lsm = sm[270:]
            s_mask = np.where(lann==xx.id)
            vox_in_structure = s_mask[0].shape[0]
            _0 = np.where(lsm[s_mask]==0)[0].shape[0]
            _1 = np.where(lsm[s_mask]>=1)[0].shape[0]
            _2 = np.where(lsm[s_mask]>=2)[0].shape[0]
            _3 = np.where(lsm[s_mask]>=3)[0].shape[0]
            _4 = np.where(lsm[s_mask]>=4)[0].shape[0]
            data.append([xx.name, 'left',vox_in_structure, _0/float(vox_in_structure), _1/float(vox_in_structure), _2/float(vox_in_structure), _3/float(vox_in_structure), _4/float(vox_in_structure)])
            if xx.name == 'Crus 1': print _1
    df = pd.DataFrame(data, columns = ['structure', 'side', 'voxels_in_structure', '0_voxel', '1_voxel', '2_voxel', '3_voxel', '4_voxel'])
    df.to_csv(injdst+'_dataframe.csv')