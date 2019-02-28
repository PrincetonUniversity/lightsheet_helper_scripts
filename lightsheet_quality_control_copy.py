#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 16:20:57 2019

@author: wanglab
"""

import os, sys, shutil
print os.getcwd()
sys.path.append('/mnt/bucket/labs/wang/pisano/Python/lightsheet')
sys.path.append('/jukebox/wang/pisano/Python/lightsheet')
from tools.registration.register import elastix_command_line_call, jacobian_command_line_call, change_interpolation_order, transformix_command_line_call
from tools.utils.io import listdirfull, makedir, load_memmap_arr, load_np, listall, load_kwargs
from skimage.external import tifffile
import skimage
from tools.utils.directorydeterminer import directorydeterminer
from scipy.ndimage.interpolation import zoom
from tools.registration.transform_cell_counts import generate_transformed_cellcount, get_fullsizedims_from_kwargs, points_resample
from tools.registration.transform_list_of_points import modify_transform_files
from tools.imageprocessing.orientation import fix_contour_orientation, fix_dimension_orientation
import matplotlib.gridspec as gridspec
from tools.conv_net.functions.dilation import dilate_with_element
from skimage.morphology import ball

if __name__ == '__main__':

    #set up
    dst = '/home/wanglab/wang/pisano/tracing_output/qc/antero_no_jg'; makedir(dst)
    lst = [xx for xx in listdirfull('/home/wanglab/wang/pisano/tracing_output/antero_4x') if 'jg' not in os.path.basename(xx)]
    
    dst = '/home/wanglab/wang/pisano/tracing_output/qc/antero_only_jg'; makedir(dst)
    lst = [xx for xx in listdirfull('/home/wanglab/wang/pisano/tracing_output/antero_4x') if 'jg' in os.path.basename(xx)]
    
    cnn_transform_type = 'affine_only_reg_to_sig';#both for regwatlas, and only affine for sig adn reg #'all', 'single': don't consider reg with sig at all
    volume_transform_type = 'single';#both for regwatlas, and only affine for sig adn reg #'all', 'single': don't consider reg with sig at all
    verbose = True
    generate_registered_overlay = False #this looks bad
    generate_downsized_overlay = True #this looks better
    #fld = '/home/wanglab/wang/pisano/tracing_output/antero_4x/20170130_tp_bl6_sim_1750r_03'    
    #loop
    for fld in lst:
        try:
            print fld
            kwargs = load_kwargs(fld)
            regvol = [xx for xx in kwargs['volumes'] if xx.ch_type == 'regch'][0]
            injvol = [xx for xx in kwargs['volumes'] if xx.ch_type == 'injch'][0]
            cellvol = [xx for xx in kwargs['volumes'] if xx.ch_type == 'cellch'][0]
            dst0 = os.path.join(dst, os.path.basename(fld)); makedir(dst0)
            dst1 = os.path.join(dst0, 'elastix'); makedir(dst1)
            
            #####check cell detection (modeled from lightsheet/tools/registration/transform_cell_counts)
            #3dunet cell dataframe
            dataframe = pd.read_csv(listdirfull(os.path.join(fld, '3dunet_output/pooled_cell_measures'), '.csv')[0])
            
            #####generate a downsized version######
            if generate_downsized_overlay:
                cellvolloaded = tifffile.imread(cellvol.resampled_for_elastix_vol)
                cnn_cellvolloaded = np.zeros_like(cellvolloaded)
                zyx = dataframe[['z','y','x']].values
                #adjust for reorientation THEN rescaling, remember full size data needs dimension change releative to resample
                fullsizedimensions = get_fullsizedims_from_kwargs(kwargs) #don't get from kwargs['volumes'][0].fullsizedimensions it's bad! use this instead
                zyx = fix_contour_orientation(zyx, verbose=verbose, **kwargs) #now in orientation of resample
                zyx = points_resample(zyx, original_dims = fix_dimension_orientation(fullsizedimensions, **kwargs), resample_dims = tifffile.imread(cellvol.resampled_for_elastix_vol).shape, verbose = verbose)[:, :3]
                zyx = np.asarray([str((int(xx[0]), int(xx[1]), int(xx[2]))) for xx in load_np(zyx)])
                from collections import Counter
                zyx_cnt = Counter(zyx)
                #now overlay
                for zyx,v in zyx_cnt.iteritems():
                    z,y,x = [int(xx) for xx in zyx.replace('(','',).replace(')','').split(',')]
                    try:
                        cnn_cellvolloaded[z,y,x] = v*100
                    except Exception, e:
                        print e
                merged = np.stack([cnn_cellvolloaded, cellvolloaded, np.zeros_like(cellvolloaded)], -1)
                merged = np.swapaxes(merged, 0,2)#reorient to horizontal
                tifffile.imsave(os.path.join(dst, '{}_points_merged_resampled_for_elastix.tif'.format(os.path.basename(fld))), merged)         
            
            #EXAMPLE USING LIGHTSHEET - assumes marking centers in the 'raw' full sized cell channel. This will transform those centers into "atlas" space (in this case the moving image)
            #in this case the "inverse transform has the atlas as the moving image in the first step, and the autofluorescence channel as the moving image in the second step 
            r2s0 = [xx for xx in listall(cellvol.inverse_elastixfld, 'reg2sig_TransformParameters.0.txt') if 'cellch' in xx][0]
            r2s1 = [xx for xx in listall(cellvol.inverse_elastixfld, 'reg2sig_TransformParameters.1.txt') if 'cellch' in xx][0]
            a2r0 = [xx for xx in listall(cellvol.inverse_elastixfld, 'atlas2reg2sig/atlas2reg_TransformParameters.0.txt') if 'cellch' in xx][0]
            a2r1 = [xx for xx in listall(cellvol.inverse_elastixfld, 'atlas2reg2sig/atlas2reg_TransformParameters.1.txt') if 'cellch' in xx][0]
            if cnn_transform_type == 'all':
                transformfiles = [r2s0, r2s1, a2r0, a2r1]
            elif cnn_transform_type == 'single':
                transformfiles = [a2r0, a2r1]
            elif cnn_transform_type == 'affine_only_reg_to_sig':    
                transformfiles = [r2s0, a2r0, a2r1]
            transformfiles = modify_transform_files(transformfiles, dst = dst1)
            
            #convert points
            converted_points = generate_transformed_cellcount(dataframe, dst1, transformfiles, lightsheet_parameter_dictionary=os.path.join(fld, 'param_dict.p'), verbose=verbose)
            
            #load and convert to single voxel loc
            zyx = np.asarray([str((int(xx[0]), int(xx[1]), int(xx[2]))) for xx in load_np(converted_points)])
            from collections import Counter
            zyx_cnt = Counter(zyx)
            
            #manually call transformix..
            transformed_dst = os.path.join(dst1, 'transformed_points'); makedir(transformed_dst)
            if volume_transform_type == 'all':
                tp0 = [xx for xx in listall(os.path.dirname(cellvol.ch_to_reg_to_atlas), 'TransformParameters.0.txt') if 'sig_to_reg' in xx and 'regtoatlas' not in xx][0]
                tp1 = [xx for xx in listall(os.path.dirname(cellvol.ch_to_reg_to_atlas), 'TransformParameters.1.txt') if 'sig_to_reg' in xx and 'regtoatlas' not in xx][0]
                transformfiles = [tp0, tp1, os.path.join(fld, 'elastix/TransformParameters.0.txt'), os.path.join(fld, 'elastix/TransformParameters.1.txt')]
            elif volume_transform_type == 'single':
                transformfiles = [os.path.join(fld, 'elastix/TransformParameters.0.txt'), os.path.join(fld, 'elastix/TransformParameters.1.txt')]
            elif volume_transform_type == 'affine_only_reg_to_sig':    
                tp0 = [xx for xx in listall(os.path.dirname(cellvol.ch_to_reg_to_atlas), 'TransformParameters.0.txt') if 'sig_to_reg' in xx and 'regtoatlas' not in xx][0]
                transformfiles = [tp0, os.path.join(fld, 'elastix/TransformParameters.0.txt'), os.path.join(fld, 'elastix/TransformParameters.1.txt')]
            transformfiles = modify_transform_files(transformfiles, dst = dst1)
            transformix_command_line_call(cellvol.resampled_for_elastix_vol, transformed_dst, transformfiles[-1])
            
            #cell_registered channel
            if generate_registered_overlay:
                cell_reg = tifffile.imread(os.path.join(transformed_dst, 'result.tif'))
                cell_cnn = np.zeros_like(cell_reg)
                errors = []
                for zyx,v in zyx_cnt.iteritems():
                    z,y,x = [int(xx) for xx in zyx.replace('(','',).replace(')','').split(',')]
                    try:
                        cell_cnn[z,y,x] = v*100
                    except Exception, e:
                        print e
                        errors.append(e)
                if len(errors)>0:
                    with open(os.path.join(dst, '{}_errors.txt'.format(os.path.basename(fld))), 'a') as flll:
                        for err in errors:
                            flll.write(str(err)+'\n')
                        flll.close()
                merged = np.stack([cell_cnn, cell_reg, np.zeros_like(cell_reg)], -1)
                #reorient to horizontal
                merged = np.swapaxes(merged, 0,2)
                tifffile.imsave(os.path.join(dst, '{}_points_merged.tif'.format(os.path.basename(fld))), merged)
            
            ####CHECK REGS       
            #get regs
            print fld, '\ngenerating registration metrics, loading....'        
            #get elastix paths
            ereg = os.path.join(regvol.elastixfld, 'result.1.tif')
            einj = injvol.ch_to_reg_to_atlas
            ecell = cellvol.ch_to_reg_to_atlas
            
            #get inverse
            a2rinj = [xx for xx in listall(injvol.inverse_elastixfld, '_atlas2reg/result.1.tif') if 'injch' in xx][0]
            r2sinj = [xx for xx in listall(injvol.inverse_elastixfld, '_reg2sig/result.1.tif') if 'injch' in xx][0]
            a2rcell = [xx for xx in listall(cellvol.inverse_elastixfld, '_atlas2reg/result.1.tif') if 'cellch' in xx][0]
            r2scell = [xx for xx in listall(cellvol.inverse_elastixfld, '_reg2sig/result.1.tif') if 'cellch' in xx][0]
            
            #load
            #vols_pths = [ereg, einj, ecell, a2rinj, r2sinj, a2rcell, r2scell]
            #names = ['auto elastix', 'inj elastix', 'cell elastix', 'atlas to auto', 'auto to inj', 'atlas to auto', 'auto to cell']
            tvols_pths = [ereg, einj, ecell, a2rinj, r2sinj, r2scell]
            tnames = ['auto\nelastix', 'inj\nelastix', 'cell\nelastix', 'atlas to\nauto', 'auto to\ninj', 'auto to\ncell']
            vols_pths=[]
            names=[]
            #this ensure things don't fail...
            for i in range(len(tvols_pths)):
                if os.path.exists(tvols_pths[i]):
                    vols_pths.append(tvols_pths[i])
                    names.append(tnames[i])
                else:
                    print('\n\nMissing {}'.format(tnames[i]))
            vols = [np.swapaxes(tifffile.imread(xx), 0, 2) for xx in vols_pths]
            
            #make figure
            fig = plt.figure(figsize=(6,8))
            rows = len(names)
            cols = 5
            tick = 0
            gs1 = gridspec.GridSpec(rows,cols)
            gs1.update(wspace=0.0000025, hspace=0.000005)
            adjust_exposure = True
            for row in range(rows):
                im = np.copy(vols[row])#*100
                step = im.shape[0] / cols
                for col in range(cols):
                    ax = plt.subplot(gs1[tick])
                    tick+=1; print tick
                    #ax = plt.subplot(rows,cols,tick)        
                    immax = np.max(im[col*step:(col+1)*step], 0)
                    if adjust_exposure: immax = skimage.exposure.equalize_hist(immax, nbins=20000)            
                    ax.imshow(immax, cmap='gray')
                    ax.set_yticklabels([])
                    ax.set_xticklabels([])
                    if col == 0: ax.set_ylabel(names[row]).set_rotation(45)
            fig.suptitle(os.path.basename(fld));                
            fig.subplots_adjust(top=0.88, wspace=0, hspace=0)
            plt.savefig(os.path.join(dst, os.path.basename(fld)+'_registration_qc.pdf'), dpi=300, transparent=True)
            plt.close()
            
            #overlay of cells
            if True:
                if False:#don't reload
                    fld = '/home/wanglab/wang/pisano/tracing_output/antero_4x/20170410_tp_bl6_lob6a_ml_repro_04'
                    pth = '/home/wanglab/wang/pisano/tracing_output/qc/antero_no_jg/20170410_tp_bl6_lob6a_ml_repro_04_points_merged_resampled_for_elastix.tif'
                    kwargs = load_kwargs(fld)
                    regvol = [xx for xx in kwargs['volumes'] if xx.ch_type == 'regch'][0]
                    injvol = [xx for xx in kwargs['volumes'] if xx.ch_type == 'injch'][0]
                    cellvol = [xx for xx in kwargs['volumes'] if xx.ch_type == 'cellch'][0]
                    dst0 = os.path.join(dst, os.path.basename(fld)); makedir(dst0)
                    dst1 = os.path.join(dst0, 'elastix'); makedir(dst1)
                    tvols_pths = [ereg, einj, ecell, a2rinj, r2sinj, r2scell]
                    tnames = ['auto\nelastix', 'inj\nelastix', 'cell\nelastix', 'atlas to\nauto', 'auto to\ninj', 'auto to\ncell']
                    vols_pths=[]
                    names=[]
                    #this ensure things don't fail...
                    for i in range(len(tvols_pths)):
                        if os.path.exists(tvols_pths[i]):
                            vols_pths.append(tvols_pths[i])
                            names.append(tnames[i])
                        else:
                            print('\n\nMissing {}'.format(tnames[i]))
                    vols = [np.swapaxes(tifffile.imread(xx), 0, 2) for xx in vols_pths]
                #done this way to prevent swapping of axes
                cnn_vol = tifffile.imread(os.path.join(dst, '{}_points_merged_resampled_for_elastix.tif'.format(os.path.basename(fld))))
                vols.append(cnn_vol)
                names.append('CNN\noverlay')
                
                #make figure
                fig = plt.figure(figsize=(6,8))
                rows = len(names)
                cols = 5
                tick = 0
                gs1 = gridspec.GridSpec(rows,cols)
                gs1.update(wspace=0.0000025, hspace=0.000005)
                adjust_exposure = True
                for row in range(rows):
                    name = names[row]
                    im = np.copy(vols[row])#*100
                    step = im.shape[0] / cols
                    if name == 'CNN\noverlay':im[:,:,:,0] = dilate_with_element(im[:,:,:,0], ball(3))
                    for col in range(cols):
                        ax = plt.subplot(gs1[tick])
                        tick+=1; print tick
                        #ax = plt.subplot(rows,cols,tick)
                        if name == 'CNN\noverlay':
                            from skimage import exposure
                            immax = np.swapaxes(np.swapaxes(np.asarray([np.max(im[col*step:(col+1)*step,:,:,i], 0) for i in range(3)]),0,1),1,2)
                            bkgrn = skimage.exposure.equalize_hist(immax[:,:,1], nbins=20000)
                            ax.imshow(bkgrn,'gray', vmin=.6, vmax=(bkgrn.max()*1.0))
                            #modify colormap
                            import matplotlib as mpl
                            my_cmap = plt.cm.Reds(np.arange(plt.cm.RdBu.N))
                            my_cmap[:1,:4] = 0.0
                            my_cmap = mpl.colors.ListedColormap(my_cmap)
                            my_cmap.set_under('w')
                            points_adjusted = skimage.exposure.equalize_hist(immax[:,:,0], nbins=20000)*1000
                            ax.imshow(points_adjusted, cmap=my_cmap, alpha=0.95, vmin=points_adjusted.min(), vmax=points_adjusted.max())
                                
                        else:
                            immax = np.max(im[col*step:(col+1)*step], 0)
                            if adjust_exposure: immax = skimage.exposure.equalize_hist(immax, nbins=20000)            
                            ax.imshow(immax, cmap='gray')
                        ax.set_yticklabels([])
                        ax.set_xticklabels([])
                        if col == 0: ax.set_ylabel(name).set_rotation(45)
                fig.suptitle(os.path.basename(fld));                
                fig.subplots_adjust(top=0.88, wspace=0, hspace=0)
                plt.savefig(os.path.join(dst, os.path.basename(fld)+'_registration_and_cell_qc.pdf'), dpi=300, transparent=True)
                plt.close()
            shutil.rmtree(dst0)
        except Exception, e:
            print e
        
    #take resampled points....
    #also then make atlas overlaid with points...
    #then make into those cfos like figures....


#%%
for channel in range(img.shape[3]):  # equalizing each channel
    img[:, :, :,channel] = exposure.equalize_hist(img[:, :, :,channel])
    
    im = np.copy(cnn_vol)
    p2, p98 = np.percentile(im[:,:,:,1], (98, 99.9999))
    im[:,:,:,1] = exposure.rescale_intensity(im[:,:,:,1], in_range=(p2, p98))
    plt.imshow(np.max(im,0))
    plt.imshow(np.max(im[:,:,:,1],0))
