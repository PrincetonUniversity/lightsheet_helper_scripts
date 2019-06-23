#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:37:12 2019

@author: tpisano
"""

#transform annotation to full sized space
from tools.utils.io import listdirfull, makedir, load_memmap_arr, load_np, listall, load_kwargs
from tools.registration.register import elastix_command_line_call, jacobian_command_line_call, change_interpolation_order, transformix_command_line_call, count_structure_lister
from tools.registration.transform_list_of_points import modify_transform_files
from tools.registration.transform_cell_counts import generate_transformed_cellcount, get_fullsizedims_from_kwargs, points_resample
from scipy.ndimage.interpolation import zoom
ann = '/home/wanglab/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_100um_edge_erosion_250um_ventricular_erosion.tif'
a2r0 = [xx for xx in listall(cellvol.inverse_elastixfld) if 'atlas2reg_TransformParameters.0' in xx and 'cellch' in xx][0]
a2r1 = [xx for xx in listall(cellvol.inverse_elastixfld) if 'atlas2reg_TransformParameters.1' in xx and 'cellch' in xx][0]
r2s0 = [xx for xx in listall(cellvol.inverse_elastixfld) if 'reg2sig_TransformParameters.0' in xx and 'cellch' in xx][0]
r2s1 = [xx for xx in listall(cellvol.inverse_elastixfld) if 'reg2sig_TransformParameters.1' in xx and 'cellch' in xx][0]
aldst = os.path.join(ldst, 'ann'); makedir(aldst)
transformfiles = modify_transform_files(transformfiles=[a2r0, a2r1, r2s0, r2s1], dst=aldst)
[change_interpolation_order(xx,0) for xx in transformfiles]
transformix_command_line_call(ann, aldst, transformfiles[-1])

#now zoom out - this is heavy!
transformed_ann = os.path.join(aldst, 'result.tif')
tann = tifffile.imread(transformed_ann)
dv0,ap0,ml0 = get_fullsizedims_from_kwargs(kwargs)

ml1,ap1,dv1 = tann.shape

#scale in dv only first and rotate to hor orientation
bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)

#make memmap
annout = os.path.join(ldst, 'transformed_annotations.npy')
annout = load_memmap_arr(annout, mode='w+', dtype='float32', shape=(dv0,ap0,ml0))

#now rotate and scale each in ap and ml
for iii,zplane in enumerate(bigdvann):
    annout[iii] = zoom(zplane, (ap0/float(ap1), ml0/float(ml1)), order=0)
    if iii%50==0: 
        annout.flush()
        print iii
annout.flush()



#%%