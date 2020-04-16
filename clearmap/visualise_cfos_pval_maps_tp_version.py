
# coding: utf-8

# In[9]:


#Generate dataframe
import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, os, itertools
sns.set_style("white")
os.chdir('/home/wanglab/wang/pisano/Python/lightsheet')
from tools.utils.io import listdirfull


ndf = pd.read_csv('/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_analysis/stim_vs_ctl/voxelization_output_TP.csv')
ndf[0:20]

from skimage.external import tifffile
import SimpleITK as sitk, numpy as np, matplotlib.pyplot as plt
vol =  tifffile.imread('/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_analysis/stim_vs_ctl/pvalues.tif')
ann = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif')
atl = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif')
import os
import scipy
import scipy.stats
os.chdir('/home/wanglab/wang/pisano/Python/lightsheet')
from tools.imageprocessing.orientation import fix_orientation
from tools.utils.io import *
#pvol = fix_orientation(vol[:,:,:,1], axes=('2','0','1'))
#nvol = fix_orientation(vol[:,:,:,0], axes=('2','0','1'))
pvol = vol[:,:,:,1]
nvol = vol[:,:,:,0]
atl = fix_orientation(atl, axes=('2','0','1'))
ann = fix_orientation(ann, axes=('2','0','1'))
assert vol.shape[:-1] == ann.shape
print ann.shape
print nvol.shape

assert pvol.shape == ann.shape


#anterior is 0 index; 100-250 is good in z
from tools.registration.allen_structure_json_to_pandas import annotation_value_to_structure,  isolate_structures_return_list, consolidate_parents_structures, annotation_location_to_structure
allen_id_table = '/home/wanglab/wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx'
#ann = '/home/wanglab/wang/pisano/Python/allenatlas/annotation_25_ccf2015.nrrd'

#find all pixel IDs present in structures
#pix_values=list(np.unique(ann[100:175].ravel().astype('float64')))

#collect names of parents OR ##optionally get sublist
if False:
    #make structures to find parents
    from tools.analysis.network_analysis import make_structure_objects
    structures = make_structure_objects('/home/wanglab/wang/pisano/Python/lightsheet/supp_files/sample_cell_count_output.xlsx')
    parent_list = list(set([yy.parent[1] for xx in pix_values for yy in structures if xx == yy.idnum]))

#find counts of highest labeled areas
if False:
    olist = annotation_location_to_structure(allen_id_table, args=zip(*np.nonzero(pvol[100:175])), ann=ann[100:175])
    srt = sorted(olist, key=lambda x: x[0], reverse=True)
    parent_list = [xx[1] for xx in srt]
    #select only subset
    parent_list=parent_list[0:13]
#manually choose locations
if False: parent_list = ['Anterior cingulate area', 'Hippocampal formation', 'Infralimbic area', 'Prelimbic area', 'Orbital area',
       'Agranular insular area', 'Somatosensory areas', 'Striatum']


#generate list of structures present
#nann, lst = consolidate_parents_structures(allen_id_table, ann, parent_list, verbose=True)

#%%
def consolidate_parents_structures_cfos(id_table, ann, namelist, verbose=False, structures=False):
    '''Function that generates evenly spaced pixels values based on annotation parents

    Removes 0 from list

    Inputs:
        id_table=path to excel file generated from scripts above
        ann = allen annoation file
        namelist=list of structues names, typically parent structures*********************

    Returns:
        -----------
        nann = new array of bitdepth
        list of value+name combinations
    '''
    if type(ann) == str: ann = sitk.GetArrayFromImage(sitk.ReadImage(ann))

    df = pd.read_excel(id_table)

    #remove duplicates and null and root
    namelist = list(set(namelist))
    namelist = [xx for xx in namelist if xx != 'null' and xx != 'root']
    namelist.sort()

    #make structures to find parents
    if not structures:
        from tools.analysis.network_analysis import make_structure_objects
        structures = make_structure_objects(id_table)

    #setup
    nann = np.zeros(ann.shape).astype('uint8')
    cmap = [xx for xx in np.linspace(1,255, num=len(namelist))]

    #populate
    for i in range(len(namelist)):
        try:
            nm=namelist[i]
            s = [xx for xx in structures if xx.name==nm][0]
            if verbose: print ('{}, {} of {}, value {}'.format(nm, i, len(namelist)-1, cmap[i]))
            nann[np.where(ann==int(s.idnum))] = cmap[i]
            for ii in s.progeny:
                if ii[3] != 'null': nann[np.where(ann==int(ii[3]))] = cmap[i]
        except Exception,e:
            print nm, e
    #sitk.Show(sitk.GetImageFromArray(nann))
    #change nann to have NAN where zeros
    nann = nann.astype('float')
    nann[nann == 0] = 'nan'

    return nann, zip(cmap[:], namelist)
#%%
tab20cmap = [plt.cm.tab20(xx) for xx in range(20)]
tab20cmap_nogray = tab20cmap[:14] + tab20cmap[16:]
#sns.palplot(tab20cmap_nogray)
import matplotlib
tab20cmap_nogray = matplotlib.colors.ListedColormap(tab20cmap_nogray, name='tab20cmap_nogray')
from matplotlib import gridspec
import matplotlib.patches as mpatches
import seaborn as sns
sns.set_style('white')
#%%
from tools.analysis.network_analysis import make_structure_objects
structures = make_structure_objects(allen_id_table)
#threshold values
pvol[pvol!=0.0] = 1.0
nvol[nvol!=0.0] = 1.0
save_dst = '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/'
no_structures_to_keep = 18
positive = True
negative = False
#loop only positives
zstep = 40
colorbar_cutoff = 65#20,60 #this is related to zstep size...(it's like  apercetnage...)
rngs = range(0, 558, zstep)
for iii in range(len(rngs)-1):
    #range rng = (100,150)
    rng = (rngs[iii], rngs[iii+1])

    #positive
    if positive:
        print rng, 'positive'
        #get highest
        olist = annotation_location_to_structure(allen_id_table, zip(*np.nonzero(pvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
        srt = sorted(olist, key=lambda x: x[0], reverse=True)
        parent_list = [xx[1] for xx in srt]
        #select only subset
        parent_list=parent_list[:no_structures_to_keep]
        nann, lst = consolidate_parents_structures_cfos(allen_id_table, ann[rng[0]:rng[1]], parent_list, verbose=True, structures=structures)

        #make fig
        plt.figure(figsize=(15,18))
        ax = plt.subplot(2,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('ABA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
        ax.set_anchor('W')

        #pvals
        ax = plt.subplot(2,1,2)
        ax.set_title('Number of positively correlated voxels in AP dimension')
        #modify colormap
        import matplotlib as mpl
        my_cmap = plt.cm.viridis(np.arange(plt.cm.RdBu.N))
        #my_cmap = plt.cm.RdYlGn(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        #plot
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        #plt.imshow(np.max(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.9) #old way
        plt.imshow(np.sum(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar()
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        #plt.tight_layout()
        pdst = os.path.join(save_dst, 'positive_overlays_zstep{}'.format(zstep)); makedir(pdst)
        plt.savefig(os.path.join(pdst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()

    #negative
    if negative:
        print rng, 'negative'
        #get highest
        olist = annotation_location_to_structure(allen_id_table, zip(*np.nonzero(nvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
        srt = sorted(olist, key=lambda x: x[0], reverse=True)
        parent_list = [xx[1] for xx in srt]
        #select only subset
        parent_list=parent_list[0:no_structures_to_keep]
        nann, lst = consolidate_parents_structures_cfos(allen_id_table, ann[rng[0]:rng[1]], parent_list, verbose=True, structures=structures)

        #make fig
        plt.figure(figsize=(15,18))
        ax = plt.subplot(2,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('ABA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
        ax.set_anchor('W')

        #pvals
        ax = plt.subplot(2,1,2)
        ax.set_title('Number of negatively correlated voxels in AP dimension')
        #modify colormap
        import matplotlib as mpl
        my_cmap = plt.cm.plasma(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        #plot
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        plt.imshow(np.sum(nvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar()
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        #plt.tight_layout()
        ndst = os.path.join(save_dst, 'negative_overlays_zstep{}'.format(zstep)); makedir(ndst)
        plt.savefig(os.path.join(ndst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()

#%%
#%%make sagittal in PMA space
from skimage.external import tifffile
import SimpleITK as sitk, numpy as np, matplotlib.pyplot as plt
vol =  tifffile.imread('/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_analysis/stim_vs_ctl/pvalues.tif')
ann = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif')
atl = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif')
import os
import scipy
import scipy.stats
os.chdir('/home/wanglab/wang/pisano/Python/lightsheet')
from tools.imageprocessing.orientation import fix_orientation
from tools.utils.io import *
pma_id_table = '/home/wanglab/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx'
from tools.analysis.network_analysis import make_structure_objects
pma_structures = make_structure_objects(pma_id_table)


#%%
#register vol into pma space...
ndst = '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/pma_space'
from tools.registration.transform_list_of_points import modify_transform_files
from tools.registration.register import transformix_command_line_call, change_interpolation_order
transformfiles = ['/home/wanglab/wang/pisano/Python/atlas/tp_to_allen/aba_to_tp/TransformParameters.0.txt', '/home/wanglab/wang/pisano/Python/atlas/tp_to_allen/aba_to_tp/TransformParameters.1.txt']
transformed_dst = os.path.join(ndst, 'elastix'); makedir(transformed_dst)
transformfiles = modify_transform_files(transformfiles, transformed_dst)
[change_interpolation_order(xx, order = 0) for xx in transformfiles]
#reoriented into dv scan space
rvol = np.copy(vol)
rvol = np.swapaxes(np.swapaxes(rvol, 0,1),0,2)
#sitk.Show(sitk.GetImageFromArray(rvol))

#split up
pvol = rvol[:,:,:,1]
pvol[pvol!=0.0] = 1.0
pdst = os.path.join(transformed_dst, 'pvol'); makedir(pdst)
pfl = os.path.join(pdst, 'pvol_aba_space.tif')
tifffile.imsave(pfl, pvol)
transformix_command_line_call(src = pfl, dst=pdst, transformfile=transformfiles[-1])

nvol = rvol[:,:,:,0]
nvol[nvol!=0.0] = 1.0
ndst = os.path.join(transformed_dst, 'nvol'); makedir(ndst)
nfl = os.path.join(ndst, 'nvol_aba_space.tif')
tifffile.imsave(nfl, nvol)
transformix_command_line_call(src = nfl, dst=ndst, transformfile=transformfiles[-1])
#%%#%%PMA space sag
#keep in sag orientation
import os
import matplotlib as mpl
import scipy
import scipy.stats
os.chdir('/home/wanglab/wang/pisano/Python/lightsheet')
from tools.imageprocessing.orientation import fix_orientation
from tools.utils.io import *
pma_id_table = '/home/wanglab/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx'
from tools.analysis.network_analysis import make_structure_objects
pma_structures = make_structure_objects(pma_id_table)
ann = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_60um_erosion.tif')
atl = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif')
pvol = tifffile.imread('/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/pma_space/elastix/pvol/result.tif')
nvol = tifffile.imread('/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/pma_space/elastix/nvol/result.tif')
print ann.shape
print pvol.shape

assert pvol.shape == ann.shape
pvol[pvol!=0.0] = 1.0
nvol[nvol!=0.0] = 1.0
save_dst = '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/pma_space'
no_structures_to_keep = 18
positive = True
negative = False
#loop only positives
zstep = 40 #30
colorbar_cutoff = 65#20,60 #this is related to zstep size...(it's like  apercetnage...)
rngs = range(0, 558, zstep)
for iii in range(len(rngs)-1):
    #range rng = (100,150)
    rng = (rngs[iii], rngs[iii+1])

    #positive
    if positive:
        print rng, 'positive'
        #get highest
        olist = annotation_location_to_structure(pma_id_table, zip(*np.nonzero(pvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
        srt = sorted(olist, key=lambda x: x[0], reverse=True)
        parent_list = [xx[1] for xx in srt]
        #select only subset
        parent_list=parent_list[:no_structures_to_keep]
        nann, lst = consolidate_parents_structures_cfos(pma_id_table, ann[rng[0]:rng[1]], parent_list, verbose=True, structures=pma_structures)

        #make fig
        plt.figure(figsize=(15,18))
        ax = plt.subplot(2,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('PMA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
        ax.set_anchor('W')

        #pvals
        ax = plt.subplot(2,1,2)
        ax.set_title('Number of positively correlated voxels in AP dimension')
        my_cmap = plt.cm.viridis(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        plt.imshow(np.sum(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar()
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        pdst = os.path.join(save_dst, 'sag_positive_overlays_zstep{}'.format(zstep)); makedir(pdst)
        plt.savefig(os.path.join(pdst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()

        #make single_overlay
        alpha=0.65
        plt.figure(figsize=(15,18))
        ax = plt.subplot(1,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('PMA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=alpha, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        leg = plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        for lh in leg.legendHandles:
            r,g,b,a = lh.get_facecolor()
            lh.set_facecolor((r,g,b,alpha))
            lh.set_edgecolor((r,g,b,alpha))
        ax.set_anchor('W')
        ax = plt.subplot(1,1,1)
        ax.set_title('Number of positively correlated voxels in AP dimension')
        my_cmap = plt.cm.inferno(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        plt.imshow(np.sum(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar(orientation='horizontal', shrink=0.75)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        pdst = os.path.join(save_dst, 'sag_positive_single_overlays_zstep{}'.format(zstep)); makedir(pdst)
        plt.savefig(os.path.join(pdst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()

    #negative
    if negative:
        print rng, 'negative'
        #get highest
        olist = annotation_location_to_structure(pma_id_table, zip(*np.nonzero(nvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
        srt = sorted(olist, key=lambda x: x[0], reverse=True)
        parent_list = [xx[1] for xx in srt]
        #select only subset
        parent_list=parent_list[0:no_structures_to_keep]
        nann, lst = consolidate_parents_structures_cfos(pma_id_table, ann[rng[0]:rng[1]], parent_list, verbose=True, structures=pma_structures)

        #make fig
        plt.figure(figsize=(15,18))
        ax = plt.subplot(2,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('PMA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
        ax.set_anchor('W')

        #pvals
        ax = plt.subplot(2,1,2)
        ax.set_title('Number of negatively correlated voxels in AP dimension')
        #modify colormap
        import matplotlib as mpl
        my_cmap = plt.cm.plasma(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        #plot
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        plt.imshow(np.sum(nvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar()
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        #plt.tight_layout()
        ndst = os.path.join(save_dst, 'sag_negative_overlays_zstep{}'.format(zstep)); makedir(ndst)
        plt.savefig(os.path.join(ndst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()

#%%#%%#%%PMA space cor
import os
import scipy
import scipy.stats
os.chdir('/home/wanglab/wang/pisano/Python/lightsheet')
from tools.imageprocessing.orientation import fix_orientation
from tools.utils.io import *
pma_id_table = '/home/wanglab/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx'
from tools.analysis.network_analysis import make_structure_objects
pma_structures = make_structure_objects(pma_id_table)
ann = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_60um_erosion.tif')
atl = tifffile.imread('/home/wanglab/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif')
#cor orientation
pvol = tifffile.imread('/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/pma_space/elastix/pvol/result.tif')
nvol = tifffile.imread('/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/pma_space/elastix/nvol/result.tif')
print ann.shape
print pvol.shape

ann = fix_orientation(ann, axes=('2','0','1'))
atl = fix_orientation(atl, axes=('2','0','1'))
pvol = fix_orientation(pvol, axes=('2','0','1'))
nvol = fix_orientation(nvol, axes=('2','0','1'))



assert pvol.shape == ann.shape
pvol[pvol!=0.0] = 1.0
nvol[nvol!=0.0] = 1.0
save_dst = '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/2019_analysis/pma_space'
no_structures_to_keep = 18
positive = True
negative = False
#loop only positives
zstep = 50
colorbar_cutoff = 70#20,60 #this is related to zstep size...(it's like  apercetnage...)
rngs = range(0, 558, zstep)
for iii in range(len(rngs)-1):
    #range rng = (100,150)
    rng = (rngs[iii], rngs[iii+1])

    #positive
    if positive:
        print rng, 'positive'
        #get highest
        olist = annotation_location_to_structure(pma_id_table, zip(*np.nonzero(pvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
        srt = sorted(olist, key=lambda x: x[0], reverse=True)
        parent_list = [xx[1] for xx in srt]
        #select only subset
        parent_list=parent_list[:no_structures_to_keep]
        nann, lst = consolidate_parents_structures_cfos(pma_id_table, ann[rng[0]:rng[1]], parent_list, verbose=True, structures=pma_structures)

        #make fig
        plt.figure(figsize=(15,18))
        ax = plt.subplot(2,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('PMA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
        ax.set_anchor('W')

        #pvals
        ax = plt.subplot(2,1,2)
        ax.set_title('Number of positively correlated voxels in AP dimension')
        #modify colormap
        import matplotlib as mpl
        my_cmap = plt.cm.viridis(np.arange(plt.cm.RdBu.N))
        #my_cmap = plt.cm.RdYlGn(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        #plot
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        #plt.imshow(np.max(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.9) #old way
        plt.imshow(np.sum(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar()
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        #plt.tight_layout()
        pdst = os.path.join(save_dst, 'cor_positive_overlays_zstep{}'.format(zstep)); makedir(pdst)
        plt.savefig(os.path.join(pdst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()

        #make single_overlay
        alpha=0.65
        plt.figure(figsize=(15,18))
        ax = plt.subplot(1,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('PMA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=alpha, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        leg = plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        for lh in leg.legendHandles:
            r,g,b,a = lh.get_facecolor()
            lh.set_facecolor((r,g,b,alpha))
            lh.set_edgecolor((r,g,b,alpha))
        ax.set_anchor('W')
        ax = plt.subplot(1,1,1)
        ax.set_title('Number of positively correlated voxels in AP dimension')
        my_cmap = plt.cm.inferno(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        plt.imshow(np.sum(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar(orientation='horizontal', shrink=0.75)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        pdst = os.path.join(save_dst, 'cor_positive_single_overlays_zstep{}'.format(zstep)); makedir(pdst)
        plt.savefig(os.path.join(pdst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()


    #negative
    if negative:
        print rng, 'negative'
        #get highest
        olist = annotation_location_to_structure(pma_id_table, zip(*np.nonzero(nvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
        srt = sorted(olist, key=lambda x: x[0], reverse=True)
        parent_list = [xx[1] for xx in srt]
        #select only subset
        parent_list=parent_list[0:no_structures_to_keep]
        nann, lst = consolidate_parents_structures_cfos(pma_id_table, ann[rng[0]:rng[1]], parent_list, verbose=True, structures=pma_structures)

        #make fig
        plt.figure(figsize=(15,18))
        ax = plt.subplot(2,1,1)
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_title('PMA structures')
        mode = scipy.stats.mode(nann, axis=0, nan_policy='omit') #### THIS IS REALLY IMPORTANT
        most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
        ann_mode = mode[0][0]
        masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
        im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
        patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
        ax.set_anchor('W')

        #pvals
        ax = plt.subplot(2,1,2)
        ax.set_title('Number of negatively correlated voxels in AP dimension')
        #modify colormap
        import matplotlib as mpl
        my_cmap = plt.cm.plasma(np.arange(plt.cm.RdBu.N))
        my_cmap[:colorbar_cutoff,:4] = 0.0
        my_cmap = mpl.colors.ListedColormap(my_cmap)
        my_cmap.set_under('w')
        #plot
        fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap='gray', alpha=1)
        plt.imshow(np.sum(nvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
        plt.colorbar()
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        ax.set_anchor('W')
        #plt.tight_layout()
        ndst = os.path.join(save_dst, 'cor_negative_overlays_zstep{}'.format(zstep)); makedir(ndst)
        plt.savefig(os.path.join(ndst, 'cfos_z{}-{}.pdf'.format(rng[0],rng[1])), dpi=300, transparent=True)
        plt.close()
#%%



#%%STOPPED
#%%
#%%
# # second coronal chunk (175 - 250)

# In[ ]:


# ## QC

# In[ ]:


def norm(im):
    up = np.percentile(im, 99.9)
    down = np.percentile(im, 0.1)
    im = np.copy(im)
    im[im < down] = down
    im[im>up] = up
    high = np.max(im)
    low = np.min(im)

    return (im - low) / (high - low)


save = '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/qc/'
lst = ['/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk01', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk02', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk03', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk04', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk05', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk06', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk07', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk08', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk10', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_mk11', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp01', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp02', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp05', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp06', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp07', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp08', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tp09', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tpal', '/home/wanglab/wang/pisano/tracing_output/cfos/201701_cfos/clearmap_analysis/bkgd5_cell105_v2_analysis/clearmap_cluster_files/201701_tpbe']

for xx in lst:
    print xx
    #im = norm(np.swapaxes(tifffile.imread(xx+'/clearmap_cluster_output/cfos_resampled.tif'), 0, 2))
    im = np.swapaxes(tifffile.imread(xx+'/clearmap_cluster_output/cfos_resampled.tif'), 0, 2)*100
    plt.figure(figsize=(15,5))
    for i in range(5):
        step = im.shape[0] / 5
        ax = plt.subplot(1,5,i+1)
        plt.imshow(np.max(im[i*step:(i+1)*step], axis=0), cmap='gray')
        ax.axis('off');
    plt.title(xx[xx.rfind('/')+1:]);
    plt.savefig(save+xx[xx.rfind('/')+1:], dpi=300)
    plt.close()
    #raw_input('Press enter to continue')
