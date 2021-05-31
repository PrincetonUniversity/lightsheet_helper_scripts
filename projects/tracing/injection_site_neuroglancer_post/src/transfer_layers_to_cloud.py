import os, json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from PIL import Image
import tifffile

from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch

from taskqueue import LocalTaskQueue
import igneous.task_creation as tc
import time

from google.cloud import storage

os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = '/home/ahoag/.cloudvolume/secrets/wanglab-pma-google-secret.json'

def upload_layer(local_layer_path,gs_layer_path):

	local_info_path = os.path.join(local_layer_path,'info')
	gs_info_path = os.path.join(gs_layer_path,'info')
	blob = bucket.blob(gs_info_path)
	
	# First make sure we haven't already uploaded this layer
	blob_exists = storage.Blob(bucket=bucket, name=gs_info_path).exists(storage_client)
	
	if blob_exists:
		print(f"Already uploaded this layer: {gs_info_path}")
		return
	blob.upload_from_filename(local_info_path)
	print(f"Uploaded local info file to cloud bucket: {gs_info_path}")
	time.sleep(2) 
	
	# Copy local cloudvolume to gs cloudvolume -> This transfers precomputed files
	localvol = CloudVolume(f'file://{local_layer_path}',parallel=True)	
	cloudvol = CloudVolume(f'gs://{bucket_name}/{gs_layer_path}',parallel=True)
	image = localvol[...]
	cloudvol[...] = image
	
	# Copy mesh
	local_mesh_path = os.path.join(local_layer_path,"mesh_mip_0_err_40")
	for basename in os.listdir(local_mesh_path):
		# make a blob (cloud file) for each localfile
		local_path = os.path.join(local_mesh_path,basename)
		gs_path = os.path.join(gs_layer_path,'mesh_mip_0_err_40',basename)
		blob = bucket.blob(gs_path)
		blob.upload_from_filename(local_path)
	return

if __name__ == '__main__':
	storage_client = storage.Client.from_service_account_json(
		'/home/ahoag/.cloudvolume/secrets/wanglab-pma-google-secret.json')
	bucket_name = 'wanglab-pma'
	bucket = storage_client.get_bucket(bucket_name)

	home_dir = '/home/ahoag/progs/pisano_etal_injections'
	viz_dir = os.path.join(home_dir,'precomputed')
	datasets = ['HSV-H129_Disynaptic','HSV-H129_Trisynaptic','PRV_Disynaptic']

	for dataset in datasets:
		print(dataset)
		dataset_dir = os.path.join(viz_dir,dataset)
		# # First transfer the merged layer for this dataset
		# merged_layer_name = f'{dataset}_merged'
		# merged_layer_dir = os.path.join(viz_dir,dataset,merged_layer_name)
		# merged_gs_layer_dir = f'injection/{dataset}/{dataset}_merged'
		# upload_layer(local_layer_path=merged_layer_dir,gs_layer_path=merged_gs_layer_dir)
		
		# # Now transfer the heatmap layer for this dataset
		# heatmap_layer_name = f'{dataset}_heatmap'
		# heatmap_layer_dir = os.path.join(viz_dir,dataset,heatmap_layer_name)
		# heatmap_gs_layer_dir = f'injection/{dataset}/{dataset}_heatmap'
		# upload_layer(local_layer_path=heatmap_layer_dir,gs_layer_path=heatmap_gs_layer_dir)



		# Next transfer all individual layers
		layer_names = sorted([x for x in os.listdir(dataset_dir) if x!='progress_dirs' and 'cells' not in x])
		for layer_name in layer_names:
			print(layer_name)
			layer_dir = os.path.join(dataset_dir,layer_name)
			assert os.path.exists(layer_dir)
			gs_layer_path = f'injection/{dataset}/{layer_name}'
			upload_layer(local_layer_path=layer_dir,gs_layer_path=gs_layer_path)
			
			
