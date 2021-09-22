#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:04:02 2020

@author: wanglab
"""

import os, sys, tifffile
import matplotlib.pyplot as plt

if __name__ == "__main__":

	#takes 1 command line args

	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	project_dir =  os.path.join(output_rootpath,
		request_name,sample_name,imaging_request,
		"rawdata/resolution_3.6x")
	
	diagnostic_plot_dir = os.path.join(project_dir,"diagnostic_plots")
	savedir = os.path.join(diagnostic_plot_dir,"downsized")
	os.makedirs(savedir,exist_ok=True)

	downsized_ch488_file = os.path.join(project_dir,"Ex_488_Em_0_downsized","downsized_for_atlas_ch488.tif")
	downsized_ch642_file = os.path.join(project_dir,"Ex_642_Em_2_downsized","downsized_for_atlas_ch642.tif")

	# Make sure planes are sagittal and look good
	channels=["488","642"]
	for channel in channels:
		downsized_file = eval(f"downsized_ch{channel}_file")
		downsized_vol = tifffile.imread(downsized_file)
		n_z_planes = downsized_vol.shape[-1] # these should be sagittal now
		z_planes = list(range(0,n_z_planes,100))
		# z_planes=[400]
		for z in z_planes:
			zstr = f'Z{str(z).zfill(4)}'
			im = downsized_vol[z]
			mean = im.mean()
			std = im.std()
			fig = plt.figure()
			ax=fig.add_subplot(1,1,1)
			ax.imshow(im,vmin=0,vmax=mean+std*3)
			ax.set_title(f'Channel {channel}, downsized sagittal plane z={z}')
			savename = os.path.join(savedir,f"channel_{channel}_{zstr}_downsized.png")
			plt.savefig(savename,format='png')
			print(f"Saved {savename}") 
			plt.clf()