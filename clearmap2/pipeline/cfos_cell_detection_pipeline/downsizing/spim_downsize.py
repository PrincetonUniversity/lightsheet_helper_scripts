#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:04:02 2020

@author: wanglab
"""

import os, numpy as np, tifffile as tif, SimpleITK as sitk, cv2, multiprocessing as mp, shutil, sys
from scipy.ndimage import zoom

cwd = os.getcwd()
utils_fld = os.path.join(cwd,"utils")
sys.path.append(utils_fld)
from pipeline_utils import fast_scandir

def resize_helper(img, dst, resizef):
	print(os.path.basename(img))
	savename = os.path.join(dst, os.path.basename(img))
	if os.path.exists(savename):
		return
	im = sitk.GetArrayFromImage(sitk.ReadImage(img))
	y,x = im.shape
	yr = int(y/resizef); xr = int(x/resizef)
	im = cv2.resize(im, (xr, yr), interpolation=cv2.INTER_LINEAR)
	tif.imsave(savename,im.astype("uint16"), compress=1)
	return "success"

if __name__ == "__main__":
	#takes 1 command line args
	output_rootpath = '/jukebox/wang/ahoag/for_cz/clearmap2_test_output'

	sample_dir = sys.argv[1].strip().rstrip("/")
	n_cores = os.cpu_count()
	# Princeton Mouse Atlas
	atlpth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"

	request_name,sample_name = sample_dir.split('/')[-2:]
	src_dir = os.path.join(sample_dir,
		'imaging_request_1/rawdata/resolution_3.6x')
	dst_dir = os.path.join(output_rootpath,request_name,
		sample_name,'imaging_request_1/rawdata/resolution_3.6x')
	channels = ['488','642']
	smartspim_prefix_dict = {
		'488':"Ex_488_Em_0",
		'642':"Ex_642_Em_2",
		}
	array_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
	channel = channels[array_id]
	smartspim_prefix = smartspim_prefix_dict[channel]
	print(f"Downsizing Channel: {channel}")
	corrected_dir_toplevel = os.path.join(src_dir,f"{smartspim_prefix}_corrected")
	try:
		corrected_dir = fast_scandir(corrected_dir_toplevel)[-1] #gets to the end of directory tree
	except:
		corrected_dir = corrected_dir_toplevel #if images already in first directory
	print("\nPath to stitched images: %s\n" % corrected_dir)
	#path to store downsized images
	dst_dir_downsized = os.path.join(dst_dir,f"{smartspim_prefix}_downsized")
	print("\nPath to storage directory: %s\n\n" % dst_dir_downsized)
	
	if not os.path.exists(dst_dir_downsized):
		os.mkdir(dst_dir_downsized)
	dst_dir_downsized_planes = os.path.join(dst_dir_downsized,"downsized_planes")
	if not os.path.exists(dst_dir_downsized_planes):
		os.mkdir(dst_dir_downsized_planes)
		
	imgs = [os.path.join(corrected_dir, xx) for xx in os.listdir(corrected_dir) if xx.endswith("tif") ]
	resizef = 5 #factor to downsize imgs by
	iterlst = [(img, dst_dir_downsized_planes, resizef) for img in imgs]
	
	p = mp.Pool(n_cores)
	p.starmap(resize_helper, iterlst)
	p.terminate()
	
	#now downsample to 140% of atlas in x,y,z
	imgs_downsized = [os.path.join(dst_dir_downsized, xx) for xx in os.listdir(dst_dir_downsized) if "tif" in xx]
	imgs_downsized.sort()
	z = len(imgs_downsized)
	y,x = sitk.GetArrayFromImage(sitk.ReadImage(imgs_downsized[0])).shape
	arr = np.zeros((z,y,x))
	
	atl = sitk.GetArrayFromImage(sitk.ReadImage(atlpth))
	atlz,atly,atlx = atl.shape #get shape, sagittal
	#read all the downsized images into a single volume 
	for i,img in enumerate(imgs_downsized):
		if i%5000==0: print(i)
		arr[i,:,:] = sitk.GetArrayFromImage(sitk.ReadImage(img)) #horizontal
	#switch to sagittal
	arrsag = np.swapaxes(arr,2,0)
	z,y,x = arrsag.shape
	print((z,y,x))
	print("\n**********downsizing....heavy!**********\n")
	
	arrsagd = zoom(arrsag, ((atlz*1.4/z),(atly*1.4/y),(atlx*1.4/x)), order=1)
	# shutil.rmtree(dst)
	# os.mkdir(dst)
	tif.imsave(os.path.join(dst_dir_downsized, f"downsized_for_atlas_ch{channel}.tif"), arrsagd.astype("uint16"))