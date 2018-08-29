#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 11:08:42 2017

@author: tpisano
"""

#to register aba to ls
import shutil, numpy as np, cv2
import subprocess as sp
import pandas as pd
from tools.registration.register import elastix_command_line_call, transformix_command_line_call
from tools.utils.io import convert_to_mhd, listdirfull, makedir #<-not using now but consider using

ls = '/jukebox/wang/pisano/Python/atlas/sagittal_atlas_20um_iso.tif'
ls_ann = '/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif'
aba = '/jukebox/wang/pisano/Python/allenatlas/average_template_25_sagittal_forDVscans.tif'
aba_ann = '/jukebox/wang/pisano/Python/allenatlas/annotation_25_ccf2015_forDVscans.tif'
out = '/jukebox/wang/pisano/Python/atlas/alignment/ls_to_aba'; makedir(out)
p0 = '/jukebox/wang/pisano/Python/atlas/alignment/parameterfolder/Order1_Par0000affine.txt'
p1 = '/jukebox/wang/pisano/Python/atlas/alignment/parameterfolder/Order2_Par0000bspline.txt'
#%%
elastix_command_line_call(fx = ls, mv=aba, out=out, parameters=[p0,p1])

#%%APPLYING TRANFORM TO GENERATE LS ANN
#need to change interpolation to zero
from tools.registration.register import change_interpolation_order
from skimage.external import tifffile
nout = '/jukebox/wang/pisano/Python/atlas/alignment/ann_transform'
makedir(nout)

#moved in new transform into nout. Then changed the output format from 16bit tiff to float32
#shutil.copy('/jukebox/wang/pisano/Python/atlas/alignment/ls_to_aba/TransformParameters.0.txt', '/jukebox/wang/pisano/Python/atlas/alignment/ann_transform/TransformParameters.0.txt')
#shutil.copy('/jukebox/wang/pisano/Python/atlas/alignment/ls_to_aba/TransformParameters.1.txt', '/jukebox/wang/pisano/Python/atlas/alignment/ann_transform/TransformParameters.1.txt')

#MANUALLY CHANGED TP1:
    #FROM (InitialTransformParametersFileName "/jukebox/wang/pisano/Python/atlas/alignment/ls_to_aba/TransformParameters.0.txt")
    #TO (InitialTransformParametersFileName "/jukebox/wang/pisano/Python/atlas/alignment/ann_transform/TransformParameters.0.txt")
    #ADDED AT END: (ResultImagePixelType "float") REPLACING: ResultImagePixelType 'SHORT'

#THIS IS REALLY IMPORTANT
#round 2 of 0 vs 1 interpolation****
[change_interpolation_order(xx, 0) for xx in listdirfull(nout) if 'TransformParameters.' in xx]
tf = [xx for xx in listdirfull(nout) if 'TransformParameters.' in xx]; tf.sort()
transformix_command_line_call(src = aba_ann, dst=nout, transformfile = tf[-1])
#change back to 3, just to keep original
#tf = [change_interpolation_order(xx, 3) for xx in listdirfull(out) if 'TransformParameters.' in xx]; tf.sort()
ls_ann = '/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif'
tifffile.imsave(ls_ann, tifffile.imread('/jukebox/wang/pisano/Python/atlas/alignment/ann_transform/result.tif'))

#%% ABA now determine volumes

#check overlay and add to it - plus delinate cri from crii
from tools.registration.transform import points_resample, points_transform, transformed_pnts_to_allen_helper_func, count_structure_lister
from skimage.external import tifffile
import collections

#load
aba_ann = tifffile.imread(aba_ann)
z,y,x = np.nonzero(aba_ann)
vals = aba_ann[z,y,x]

#determine each px value and convert to structures
xlfl = '/jukebox/wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx'
counts = count_structure_lister(xlfl, *vals)

#adjust columns voxel counts and save out (check if cell_count or count)
df = counts.rename(columns={'cell_count': 'voxels_in_structure'})
#df.insert(8, 'count', 0)
xlout = '/jukebox/wang/pisano/Python/lightsheet/supp_files/allen_id_table_w_voxelcounts.xlsx'
df.to_excel(xlout)

#modify cell count:
cc = '/jukebox/wang/pisano/Python/lightsheet/supp_files/sample_cell_count_output.xlsx'
cdf = pd.read_excel(cc)
cdf.insert(7, 'voxels_in_structure', df['voxels_in_structure'])
cdf.to_excel(cc)

#%%
#%% LS now determine volumes
#check overlay and add to it - plus delinate cri from crii
from tools.registration.transform import points_resample, points_transform, transformed_pnts_to_allen_helper_func, count_structure_lister
from skimage.external import tifffile
import collections

#load
ls_ann = tifffile.imread(ls_ann)
z,y,x = np.nonzero(ls_ann)
vals = ls_ann[z,y,x]

#gut check to make sure all pixels map!
#xx = transformed_pnts_to_allen_helper_func(np.asarray(zip(z,y,x)).reshape(len(z),3),ls_ann, xlfl)

#determine each px value and convert to structures
xlfl = '/jukebox/wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx'
counts = count_structure_lister(xlfl, *vals)

#adjust columns voxel counts and save out (check if cell_count or count)
df = counts.rename(columns={'cell_count': 'voxels_in_structure'})
#df.insert(8, 'count', 0)
xlout = '/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx'
df.to_excel(xlout)


#%%
#calculate delta in volume / per structure - this gives another "goodness of registration"
ls = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx') 
aba = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/allen_id_table_w_voxelcounts.xlsx')
aba_scale = 25
ls_scale = 20

fig = plt.figure()
ax = plt.subplot(111)
ax.set_yscale=('log')
a=[]; l=[]; d = []
for i in range(len(aba)):
    abar = aba.iloc[i]
    lsr = ls.iloc[i]
    if abar['voxels_in_structure'] > 0:
        a.append(abar['voxels_in_structure']*aba_scale)
        l.append(lsr['voxels_in_structure']*ls_scale)
        d.append(abs(abar['voxels_in_structure']*aba_scale - lsr['voxels_in_structure']*ls_scale))
    
plt.scatter([0]*len(a),a, c='r')
plt.scatter([1]*len(l),l, c='b')
plt.xticks((0,1), ['ABA', 'LS'])


#change in volume of ~5%:
deltavol = abs(np.sum(a)-np.sum(l))/ float(np.sum(l))

#look to see if removing the cb does this account for it...
#get posterior pix ids
tann = ls_ann[:,500:,:]
z,y,x = np.nonzero(tann)
pixids = collections.Counter(tann[z,y,x]).values()

a=[]; l=[]; d = []
for i in range(len(aba)):
    abar = aba.iloc[i]
    lsr = ls.iloc[i]
    if abar['voxels_in_structure'] > 0 and lsr.id not in pixids:
        a.append(abar['voxels_in_structure']*aba_scale)
        l.append(lsr['voxels_in_structure']*ls_scale)
        d.append(abs(abar['voxels_in_structure']*aba_scale - lsr['voxels_in_structure']*ls_scale))
deltavol = abs(np.sum(a)-np.sum(l))/ float(np.sum(l))

#^some logic is messed up above?



#%%
data = np.random.random([10, 3, 4])
for idx,val in np.ndenumerate(data): # generator, so can be used only once
    z,y,x = idx


#%%
#################################################################
#WORKING ON FINISHING UP MAKING THE ATLAS...SEEMS THAT THERE ARE UNANNOTATED THINGS. NEED TO CHECK THIS USING SCRIPT BELOW

####NEED TO NOW GET DECIMALS TO ROUND UP OR DOWN....

#load
aba_ann = tifffile.imread('/jukebox/wang/pisano/Python/allenatlas/annotation_25_ccf2015_forDVscans.tif')#.astype('int32')
tp_ann = tifffile.imread('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif')#.astype('int32')
amira_annotations = tifffile.imread('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso.crura+vermlabels.tif')

#dtypes:
for i in [aba_ann, tp_ann, amira_annotations]:
    print i.dtype
    
df = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx')
lst = df.atlas_id.tolist()
lst = df.id.tolist()

#zero out areas that are annotated and find the nonannotations
a = aba_ann.copy()
for i, atlas_id in enumerate(lst):
    a[a==atlas_id]=0
    if i%50==0: print i, len(lst)
    
unannotated = np.where(a>0)
unannotated_pixel_vals = [int(xx) for xx in np.unique(a)[1:]] #removing zero because of above, [182305696, 182305712, 312782560, 312782592, 312782656]
print len(unannotated_pixel_vals)


#test to see where unannotated_pixel_vals map?
from tools.registration.transform import count_structure_lister
df = count_structure_lister('/jukebox/wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx', *unannotated_pixel_vals)

#%%
#merge tp and amira
aba_ann = tifffile.imread('/jukebox/wang/pisano/Python/allenatlas/annotation_25_ccf2015_forDVscans.tif')
df = pd.read_excel('/home/wanglab/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')

'''
#########################
#####IF RUNNING FROM SCRATCH NEED THIS TO PUT IN THE NEW ROWS FOR DECLIVE#####
df = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')
ovly = tifffile.imread('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif')
lob6 = df[df.name.str.contains('Declive')]
cols = df.columns
#set new lines, ****DO THIS ONLY ONCE!!!**
if len(df[df.name=='Declive (VI), subdivision A']) == 0:
    line0 = pd.DataFrame(data = [['Declive (VI)', 'DEC', 1138, '399', '645', 'Declive (VI)', 'DEC', 0]], columns = cols, index=[lob6.index[0]+1])
    line1 = pd.DataFrame(data = [['Declive (VI), subdivision A', 'DEC6a', 936, '399a', '1138', 'Declive (VI)', 'DEC', len(np.where(ovly==936)[0])]], columns = cols, index=[lob6.index[0]+1])
    line2 = pd.DataFrame(data = [['Declive (VI), subdivision B', 'DEC6b', new_mapping['l6b'], '399b', '1138', 'Declive (VI)', 'DEC', len(np.where(ovly==new_mapping['l6b'])[0])]], columns = cols, index=[lob6.index[0]+1])
    #combine
    df2 = pd.concat([df.ix[:lob6.index[0]-1], line0, line1, line2, df.ix[lob6.index[0]+1:]]).reset_index(drop=True)
    df2.to_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')
#########################    
'''


#aba_ann = np.round(aba_ann)
tp_ann = tifffile.imread('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif')#.astype('int32')
#tp_ann = np.round(tp_ann)
#amira_annotations = tifffile.imread('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso.crura+vermlabels.tif')
#issues with above so opened it in imagej, converted to rgb color, split, merged, and saved as 8bit
amira_annotations = tifffile.imread('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso.crura+vermlabels_rgbcolor_to_8bit.tif')
new = tp_ann.copy()

#find unique ~900 ids for use. USE ID NOT ATLAS_ID
def missing_elements(L):
    start, end = L[0], L[-1]
    return sorted(set(range(start, end + 1)).difference(L))
ndf = df[df.id.isin(range(700,1250))]
ununsed_ids = missing_elements(list(set(ndf.id.unique()))) #[1013, 1115, 1122, 1130, 1134, 1135, 1136, 1137, 1138]
#Left = 0-270, right = 270-540
new_mapping = {xx:yy for xx,yy in zip(['left_cri', 'left_crii', 'right_cri', 'right_crii', 'l6b', 'l7', 'l8', 'l9'], [1056, 1064, 1056, 1064, 1134, 944, 951, 957])} #944=l7, 

#amira_annotations = amira_annotations.astype('uint8')
amira_annotations = np.swapaxes(np.swapaxes(amira_annotations, 1,2), 2,3)
###MAPPING IS OFF, NEED TO FIX HERE***** change to uint16 ^ and learn mapping (this is because ImageJ can't import uint16 color...)

#mapping, l6a will be all remaining 6
dct = {'left_cri' : (0,0,255),
       'left_crii' :  (255,255,0),
       'right_cri' : (255,0,0),
       'right_crii' : (0,255,0),
       'l6b' : (127,0,127),
       'l7' : (0,127,127),
       'l8' : (0,255,255), 
       'l9' : (255,0,255)}

#find pixel id locs
for area, pixid in dct.iteritems():
    #zyx = zip(*np.where(amira_annotations==pixid)[:3]); #dct[area] = {pixid: } #zyx
    #for some reason the 127 doesn't work?
    zyx = np.where((amira_annotations[:,:,:,0]==pixid[0])&(amira_annotations[:,:,:,1]==pixid[1])&(amira_annotations[:,:,:,2]==pixid[2]))[:3]
    print area, pixid, new_mapping[area], zyx[0].shape
    new[zyx] = new_mapping[area]
#save
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif', new)    
import cv2
#look for differences in l6, cri etc...
ovly = tifffile.imread('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif')
df = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')
#dfs
lob6 = df[df.name.str.contains('Declive')]
un6 = np.where(ovly==lob6.id.values[0])

#test loading fiji rois
from tools.conv_net.input.read_roi import read_roi_zip
from tools.registration.transform import swap_cols
#roi notes: annotating lob6B**** and YOU MUST DO and 'AND' with LOB6 voxels, annotated by opening up sag atlas and sag lob6.tif above
roipth = '/jukebox/wang/pisano/Python/atlas/amira/lob6_annotations/lob6bRoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]==df[df.name.str.contains('Declive')].id.values[0]) & (blank>0))
    ovly[z, yy, xx] = new_mapping['l6b']
    
#make everything remaining 6a
lob6 = df[df.name.str.contains('Declive')]
un6 = np.where(ovly==lob6.id.values[0])
#not done yet, just need to update df
#ovly[un6] = new_mapping['l6a']

#do this for all to check!
tmp = np.zeros_like(ovly).astype('uint8')
tmp[ovly==lob6.id.values[0]]=255
tmp[ovly==new_mapping['l6b']]=175
tmp[ovly==new_mapping['l7']]=150
tmp[ovly==new_mapping['l8']]=125
tmp[ovly==new_mapping['l9']]=100
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/amira/lob6.tif',tmp) 
#
'''don't need l7 or l8 because I just added in new ones, don't need to parse by 6a/b or r vs l
lob7 = df[df.name.str.contains('Folium-tuber vermis')]
print un7[0].shape
un7 = np.where(ovly==lob7.id.values[0])
lob8 = df[df.name.str.contains('Pyramus')]
un8 = np.where(ovly==lob8.id.values[0])
print un8[0].shape
lob9 = df[df.name.str.contains('Uvula')]
'''
#find voxels that are Left = 0-270, right = 270-540
ansi = df[df.name.str.contains('Ansiform')] #find where I missed
cri = df[df.name.str.contains('Crus 1')] #find where I missed
crii = df[df.name.str.contains('Crus 2')]
#zyx = tuple([np.concatenate((xx,yy,zz),axis=0) for xx,yy,zz in zip(np.where(ovly==crura.id.values[0]), np.where(ovly==crura.id.values[1]), np.where(ovly==crura.id.values[5]))])
zyx = np.where(ovly==ansi.id.values[0])
print len(zyx[0])



tifffile.imsave('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif', ovly)   
ovly = tifffile.imread('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif')
#test loading fiji rois
from tools.conv_net.input.read_roi import read_roi_zip
from tools.registration.transform import swap_cols
#roi notes: annotating lob6B**** and YOU MUST DO and 'AND' with LOB6 voxels, annotated by opening up sag atlas and sag lob6.tif above
roipth = '/jukebox/wang/pisano/Python/atlas/amira/crus_annotations/cri_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]==df[df.name.str.contains('Ansiform')].id.values[0]) & (blank>0))
    ovly[z, yy, xx] = new_mapping['left_cri']
    

#set PM annotations
roipth = '/jukebox/wang/pisano/Python/atlas/amira/pm_annotations/pm_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]==df[df.name.str.contains('Ansiform')].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name.str.contains('Paramedian lobule')].id.values[0]
    
#fix simplex being overlwhemed by l6a annotations
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/l6_overwhelming_simplex_annotations/l6_overwhelm_simRoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]==df[df.name=='Declive (VI), subdivision A'].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name.str.contains('Simple')].id.values[0]
    
#minor trim now simplex replacing w 6a
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/trim_l6a_w_simplex/triml6_w_simplexRoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name.str.contains('Simple')].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name=='Declive (VI), subdivision A'].id.values[0]    


#minor trim l6a from crI
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/remove_6a_from_cri/rm6a_from_cri_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name=='Declive (VI), subdivision A'].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name.str.contains('Crus 1')].id.values[0]    
    
#fix lb6b removing some of 6a
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/6b_fix/6b_remove6a_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name=='Declive (VI), subdivision A'].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = new_mapping['l6b']
    
#fix lb7 removing some of 6a
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/7fix/7_remove6a_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name=='Declive (VI), subdivision A'].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = new_mapping['l7']


#fix pm removing some of crii
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/pm_fix_from_crii/pm_rmcrII_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name.str.contains('Crus 2')].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name.str.contains('Paramedian lobule')].id.values[0]

#fix pm removing some of crii
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/7_from_crii/7_from_crii_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name.str.contains('Crus 2')].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = new_mapping['l7']

#fix cp removing some of pm
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/cop_pyramids_from_pm/cop_py_from_pm_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name.str.contains('Paramedian lobule')].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name.str.contains('Copula pyramidis')].id.values[0]
    
#fix cii removing some of 6a/b
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/crii_from_6a/crii_from6ab_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where(((ovly[z]== df[df.name=='Declive (VI), subdivision A'].id.values[0])|(ovly[z]== df[df.name=='Declive (VI), subdivision B'].id.values[0])) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name.str.contains('Crus 2')].id.values[0]

#fix 9 removing some of cp
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/9_from_cp/9_from_cp_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where(((ovly[z]== df[df.name=='Copula pyramidis'].id.values[0])|(ovly[z]== df[df.name=='Copula pyramidis'].id.values[0])) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name.str.contains('Uvula')].id.values[0]

#finally some cp removing some of 9
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/cp_from_9/cp_from_9_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name.str.contains('Uvula')].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name=='Copula pyramidis'].id.values[0]


#finally some cp removing some of 9
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/cp_from_8/cp_from_8_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== df[df.name.str.contains('Pyramus')].id.values[0]) & (blank>0)) #it is correct in being ansiform
    ovly[z, yy, xx] = df[df.name=='Copula pyramidis'].id.values[0]



#finally remaining set as crii
ovly[ovly == ansi.id.values[0]] = new_mapping['left_crii']


#optional inspect
tmp = np.zeros_like(ovly).astype('uint8')
tmp[ovly==ansi.id.values[0]]=250
tmp[ovly==cri.id.values[0]]=200
tmp[ovly==crii.id.values[0]]=175
tmp[ovly==df[df.name.str.contains('Paramedian lobule')].id.values[0]]=150
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/amira/crura.tif',tmp) 
    
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif', ovly)    
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif', ovly) #<--final atlas   


#find all 0 pixels on annotation that are nonzero on atlas. Find the closest nonzero pixel and set as that. Use dist func
ovly = tifffile.imread('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif')
ovlyc = np.copy(ovly)
atl = tifffile.imread('/home/wanglab/wang/pisano/Python/atlas/sagittal_atlas_20um_iso.tif')
atl[atl<=3]=0
tifffile.imsave('/home/wanglab/wang/pisano/Python/atlas/sagittal_atlas_20um_iso.tif', atl)
import scipy
#find all nonzero annotations and put into list
ovly_nz = np.swapaxes(np.asarray(np.where(ovly>0)),0,1)
#find all nonaccounted for pixels (pixels that are nonzero in atlas but don't have annotations)
z,y,x=np.where((ovly==0)&(atl>0))
a = scipy.spatial.KDTree(np.asarray(ovly_nz))
for pnt in np.asarray(zip(z,y,x)):
    #find the closest point
    dist, idx = a.query(pnt)
    #find the annotation value of that point and insert it into the orignal annotation
    ovly[tuple(pnt)] = ovlyc[tuple(ovly_nz[idx])]
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif', ovly)   

#now remove annotations that don't have atlas signal
ovly = tifffile.imread('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif')
atl = tifffile.imread('/home/wanglab/wang/pisano/Python/atlas/sagittal_atlas_20um_iso.tif')
ovly[np.where((ovly>0)&(atl==0))] = 0
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif', ovly)   

#minor changes and we should be good
from tools.utils.io import change_arr_given_indices_condition
#zplane 255 (fiji, thus 254) - px value 1064 needs to become 936 and 1134
ovly = change_arr_given_indices_condition(ovly, indices = (254, slice(544, 564), slice(106, 146)), condition = 1064, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (254, slice(558, 561), slice(103, 107)), condition = 1064, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (254, slice(553, 559), slice(101,107)), condition = 1064, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (254, slice(560, 580), slice(96,112)), condition = 1064, new_value = 1134)
ovly = change_arr_given_indices_condition(ovly, indices = (254, slice(560, 580), slice(112,124)), condition = 1064, new_value = 944)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(223,240), slice(554,552), slice(149,162)), condition = 944, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(98,100), slice(553,559), slice(201,206)), condition = 1007, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(200,204), slice(550,559), slice(90,100)), condition = 1007, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(208,213), slice(550,555), slice(115,125)), condition = 1134, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (218, slice(553,555), slice(113,125)), condition = 1134, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(203,329), slice(567,615), slice(99,138)), condition = 936, new_value = 944)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(331,328), slice(546,554), slice(82,95)), condition = 1007, new_value = 936)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(220,240), slice(540,556), slice(148,162)), condition = 944, new_value = 936)


#
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/lob7_from_8/lob7_from_8.RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 951) & (blank>0))
    ovly[z, yy, xx] = 944

roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/lob8_from_7/lob8_from_7RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 944) & (blank>0))
    ovly[z, yy, xx] = 951


roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/lob9_from_8/lob(_from_8_RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 951) & (blank>0))
    ovly[z, yy, xx] = 957
    
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/7_from_6/7_from_6.RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 936) & (blank>0))
    ovly[z, yy, xx] = 944
    
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/6_from_6b_7/6_from_6b_7.RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where(((ovly[z]== 944)|(ovly[z]== 1134)) & (blank>0))
    ovly[z, yy, xx] = 936

roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/pm_from_cp/pm_from_cp.roi.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 1025) & (blank>0))
    ovly[z, yy, xx] = 1064
    
ovly = change_arr_given_indices_condition(ovly, indices = (175, slice(553,560), slice(100,110)), condition = 1007, new_value = 1056)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(197,198), slice(558,563), slice(100,110)), condition = 1134, new_value = 1056)
ovly = change_arr_given_indices_condition(ovly, indices = (320, slice(562,564), slice(319,322)), condition = 936, new_value = 1134)

roipth =  '/home/wanglab/wang/pisano/Python/atlas/amira/pm_from_cp/1.RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 1025) & (blank>0))
    ovly[z, yy, xx] = 1064
    
    
roipth =  '/home/wanglab/wang/pisano/Python/atlas/amira/cp_from_pm/cp_from_pm.RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 1064) & (blank>0))
    ovly[z, yy, xx] = 1025


roipth =  '/home/wanglab/wang/pisano/Python/atlas/amira/6b_from_6/6b_from_6.RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 936) & (blank>0))
    ovly[z, yy, xx] = 1134
    
roipth = '/home/wanglab/wang/pisano/Python/atlas/amira/lateral_fix/lateral_fix.RoiSet.zip'
rois = [xx for xx in read_roi_zip(roipth, include_roi_name=True) if '.zip.roi' not in xx[0]]
#format so each rois is [z,y,x, [contour]]; remeber IMAGEJ has one-based numerics for z plane. NOTE roi is ZYX, contour[XY???], why swap_cols
rois = [(map(int, xx[0].replace('.roi','').split('-')), xx[1]) for xx in rois]
rois = [(xx[0][0]-1, xx[0][1], xx[0][2], swap_cols(xx[1], 0,1)) for xx in rois]        
for roi in rois:
    z,y,x,yxcontour = roi
    blank = np.zeros_like(ovly[0])
    blank = cv2.fillPoly(blank, [np.int32(yxcontour)], color=255)
    yy, xx = np.where((ovly[z]== 1064) & (blank>0))
    ovly[z, yy, xx] = 944
    
    
ovly = change_arr_given_indices_condition(ovly, indices = (slice(95,106), slice(460,485), slice(116,145)), condition = 1064, new_value = 1056)
ovly = change_arr_given_indices_condition(ovly, indices = (slice(228,230), slice(520,540), slice(66,100)), condition = 1007, new_value = 936)


tifffile.imsave('/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif', ovly)   
    

'''
#DID NOT MISS ANY CRUS THEREFORE IGNORE BELOW!
lzyx = zyx[0]<270
rzyx = zyx[0]>=270
lc = tuple([xx[lzyx] for xx in zyx])
rc = tuple([xx[rzyx] for xx in zyx])
#remap
ovly[lc] = new_mapping['left_cri']
ovly[rc] = new_mapping['right_cri']
tifffile.imsave('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif', ovly)
#need to still distinguish between CrI and CrII
'''
#%%finally update df 
df = pd.read_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')
ovly = tifffile.imread('/jukebox/wang/pisano/Python/atlas/amira/sagittal_atlas_20um_iso_overlay.tif')
lob6 = df[df.name.str.contains('Declive')]
cols = df.columns
#set new lines, ****DO THIS ONLY ONCE!!!**
if len(df[df.name=='Declive (VI), subdivision A']) == 0:
    line0 = pd.DataFrame(data = [['Declive (VI)', 'DEC', 1138, '399', '645', 'Vermal regions', 'VERM', 0]], columns = cols, index=[lob6.index[0]+1])
    line1 = pd.DataFrame(data = [['Declive (VI), subdivision A', 'DEC6a', 936, '399a', '1138', 'Declive (VI)', 'DEC', len(np.where(ovly==936)[0])]], columns = cols, index=[lob6.index[0]+1])
    line2 = pd.DataFrame(data = [['Declive (VI), subdivision B', 'DEC6b', new_mapping['l6b'], '399b', '1138', 'Declive (VI)', 'DEC', len(np.where(ovly==new_mapping['l6b'])[0])]], columns = cols, index=[lob6.index[0]+1])
    #combine
    df2 = pd.concat([df.ix[:lob6.index[0]-1], line0, line1, line2, df.ix[lob6.index[0]+1:]]).reset_index(drop=True)
    df2.to_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')
#set new voxel counts for a
crura = df.ix[975:1040]
for i, row in df.iterrows():
    vox_count = np.where(ovly == row['id'])
    print row['name'], len(vox_count[0])
    df.ix[i, 'voxels_in_structure']=len(vox_count[0])
df.to_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')

#change simple to simplex on sheet
for i, row in df[df.name.str.contains('Simple')].iterrows():
    df.ix[i, 'name']=row['name'].replace('Simple ', 'Simplex')
if 'cell_count' not in df: df['cell_count'] = 0
df.to_excel('/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx')
