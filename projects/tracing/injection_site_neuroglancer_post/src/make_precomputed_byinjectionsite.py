import os, json, glob, shutil, gzip
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from PIL import Image
import tifffile

from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch

from taskqueue import LocalTaskQueue
import igneous.task_creation as tc



def make_info_file(volume_size,resolution,layer_dir,voxel_offset=[0,0,0],commit=True):
	""" 
	---PURPOSE---
	Make the cloudvolume info file.
	---INPUT---
	volume_size     [Nx,Ny,Nz] in voxels, e.g. [2160,2560,1271]
	pix_scale_nm    [size of x pix in nm,size of y pix in nm,size of z pix in nm], e.g. [5000,5000,10000]
	commit          if True, will write the info/provenance file to disk. 
					if False, just creates it in memory
	atlas_type      if provided, will add a key to the info file: 
					'atlas_type': atlas_type
	"""
	info = CloudVolume.create_new_info(
		num_channels = 1,
		layer_type = 'segmentation', # 'image' or 'segmentation'
		data_type = 'uint8', # 
		encoding = 'raw', # other options: 'jpeg', 'compressed_segmentation' (req. uint32 or uint64)
		resolution = resolution, # Size of X,Y,Z pixels in nanometers, 
		voxel_offset = voxel_offset, # values X,Y,Z values in voxels
		chunk_size = [ 1024,1024,1 ], # rechunk of image X,Y,Z in voxels -- only used for downsampling task I think
		volume_size = volume_size, # X,Y,Z size in voxels
		)

	vol = CloudVolume(f'file://{layer_dir}', info=info)
	vol.provenance.description = "Test on spock for profiling precomputed creation"
	vol.provenance.owners = ['ahoag@princeton.edu'] # list of contact email addresses
	if commit:
		vol.commit_info() # generates info json file
		vol.commit_provenance() # generates provenance json file
		print("Created CloudVolume info file: ",vol.info_cloudpath)
	return vol

def process_slice(z):
	if os.path.exists(os.path.join(progress_dir, str(z))):
		print(f"Slice {z} already processed, skipping ")
		return
	if z >= z_dim:
		print("Index {z} >= z_dim of volume, skipping")
		return
	print('Processing slice z=',z)
	array = image[z].reshape((1,y_dim,x_dim)).T
	vol[:,:, z] = array
	touch(os.path.join(progress_dir, str(z)))
	return "success"

def make_mesh(vol,cores=8):
	# Mesh on 8 cores, use parallel=True to use all cores
	cloudpath = vol.cloudpath
	with LocalTaskQueue(parallel=cores) as tq:
	  tasks = tc.create_meshing_tasks(cloudpath, mip=0, shape=(256, 256, 256))
	  tq.insert_all(tasks)
	  tasks = tc.create_mesh_manifest_tasks(cloudpath)
	  tq.insert_all(tasks)
	print("Done!")	

if __name__ == '__main__':
	home_dir = '/home/ahoag/progs/pisano_etal_injections'
	resolution = (20000,20000,20000) # 20 micron isotropic
	viz_dir = os.path.join(home_dir,'precomputed')

	data_dir = os.path.join(home_dir,'data','pooled_by_injection_site_volumes')

	datasets = ['HSV-H129_Disynaptic','HSV-H129_Trisynaptic','PRV_Disynaptic']

	for dataset in datasets:
		print(dataset)
		pooled_vol_files = glob.glob(data_dir + f'/{dataset}*tif')
		for pooled_vol_file in pooled_vol_files:
			primary_inj_site = pooled_vol_file.split(dataset+"_")[-1].split('.tif')[0]
			print(primary_inj_site)
			layer_name = f'{dataset}_pooled_by_{primary_inj_site}'			
			layer_dir = os.path.join(viz_dir,dataset,layer_name)
			if os.path.exists(layer_dir):
				print(f"{layer_dir} already exists. Skipping")
				continue
			progress_parentdir = os.path.join(viz_dir,dataset,'progress_dirs')
			merged_vol_file = os.path.join(data_dir,f'{dataset}_{primary_inj_site}.tif')
			image = tifffile.imread(merged_vol_file)
			z_dim,y_dim,x_dim = image.shape
			volume_size = (x_dim,y_dim,z_dim)		
			vol = make_info_file(	
				volume_size=volume_size,
				layer_dir=layer_dir,
				resolution=resolution,
				)
			progress_dir = mkdir(progress_parentdir + f'/progress_{layer_name}') # unlike os.mkdir doesn't crash on prexisting 

			done_files = set([ int(z) for z in os.listdir(progress_dir) ])
			all_files = set(range(vol.bounds.minpt.z, vol.bounds.maxpt.z)) 

			to_upload = [ int(z) for z in list(all_files.difference(done_files)) ]
			to_upload.sort()
			n_to_upload = len(to_upload)
			print(f"Have {n_to_upload} z planes to upload")
			with ProcessPoolExecutor(max_workers=10) as executor:
				for job in executor.map(process_slice,to_upload):
					try:
						print(job)
					except Exception as exc:
						print(f'generated an exception: {exc}')
			print("Done uploading planes")
			print()
			print("Making mesh")
			mesh_dir = os.path.join(layer_dir,'mesh_mip_0_err_40') 
			if os.path.exists(mesh_dir):
				print("Mesh already exists. Skipping")
				print()
				continue
			make_mesh(vol,cores=10)
			# Unzip the gzipped files
			print("Unzipping gzipped files")
			for gzipped_file in glob.glob(mesh_dir + '/*gz'):
				print(gzipped_file)
				unzipped_file = gzipped_file.replace('.gz','')
				with gzip.open(gzipped_file,'rb') as f_in:
					with open(unzipped_file,'wb') as f_out:
						shutil.copyfileobj(f_in,f_out)
				# now remove the original
				os.remove(gzipped_file)
			print()
