#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:01:17 2020

@author: wanglab
"""

import sys,os
sys.path.append("/jukebox/wang/ahoag/brainpipe")
from tools.registration.register import elastix_command_line_call

if __name__ == '__main__':
	#takes 6 command line arguments max
	array_id = int(os.environ["SLURM_ARRAY_TASK_ID"]) # 0 = atlas -> reg, 1 = reg -> cell

	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	src = os.path.join(output_rootpath,request_name,sample_name,
		imaging_request,"rawdata/resolution_3.6x")
	reg = "Ex_488_Em_0"
	cell = "Ex_642_Em_2"
	n_cores = os.cpu_count()
	
	param_fld = "/jukebox/wang/ahoag/brainpipe/parameterfolder" 
	atl = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif" # PMA 20 micron isotropic

	out1 = os.path.join(src, "elastix_inverse_transform")
	os.makedirs(out1, exist_ok=True)
	
	if array_id == 0:

		print("Doing inverse transform atlas -> reg")
		#atlas to registration vol
		#inverse transform
		fx1 = os.path.join(src, "Ex_488_Em_0_downsized", "downsized_for_atlas_ch488.tif")
		mv1 = atl
		assert os.path.exists(fx1)
		assert os.path.exists(mv1)
		print("\nPath to downsized vol for inverse registration to atlas: %s" % fx1)
		print("\nPath to atlas: %s" % mv1)
		
		
		params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
		#run
		e_out, transformfiles = elastix_command_line_call(fx1, mv1, out1, params)

	elif array_id == 1:
		#registration vol to cell vol
		#inverse transform
		print("Doing inverse transform reg -> cell")

		mv2 = os.path.join(src, "Ex_488_Em_0_downsized","downsized_for_atlas_ch488.tif")
		fx2 = os.path.join(src, "Ex_642_Em_2_downsized","downsized_for_atlas_ch642.tif")
		
		assert os.path.exists(fx2)
		assert os.path.exists(mv2)
		
		out2 = os.path.join(out1,"488_to_642")
		os.makedirs(out2, exist_ok=True)
					
		params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
		#run
		e_out, transformfiles = elastix_command_line_call(fx2, mv2, out2, params)
		
		
