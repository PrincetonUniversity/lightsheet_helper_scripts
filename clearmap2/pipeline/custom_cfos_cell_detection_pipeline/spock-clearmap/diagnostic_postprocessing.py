# General imports
import os,sys,tifffile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.append('/jukebox/braininit/lightsheet/ClearMap2')

# ClearMap imports
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io

if __name__ == "__main__":
	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	project_dir =  os.path.join(output_rootpath,
		request_name,sample_name,imaging_request,
		"rawdata/resolution_3.6x")
	
	diagnostic_plot_dir = os.path.join(project_dir,"diagnostic_plots")
	savedir = os.path.join(diagnostic_plot_dir,"postprocessing")
	os.makedirs(savedir,exist_ok=True)

	ws = wsp.Workspace('CellMap',directory=project_dir)
	ws.debug = False
	im_source = ws.source('stitched')
	cells_transformed_source = ws.source('cells', postfix='transformed_to_atlas')

	""" First plot: unfiltered cells: histograms of coordinates, size, intensity, etc.. """
	fig = plt.figure(); plt.clf();
	names = cells_transformed_source.dtype.names;
	nx,ny = 2,3
	for i, name in enumerate(names):
	    plt.subplot(nx, ny, i+1)
	    plt.hist(cells_transformed_source[name]);
	    plt.title(name)
	    if name == 'intensity':
	        plt.xlabel('counts')
	    elif name == 'region':
	        plt.xlabel('region_id')
	    else:
	        plt.xlabel('voxels')
	plt.tight_layout()

	savename=os.path.join(savedir,"transformed_cells_distributions_unfiltered.png")
	plt.savefig(savename,format='png')
	print(f"Saved {savename}")
	sys.stdout.flush()

	""" Next plot: filtered cells: histograms of coordinates, size, intensity, etc.. """
	cells_transformed_filtered_source = ws.source('cells', postfix='transformed_to_atlas_filtered_20px')

	fig = plt.figure(); plt.clf();
	names = cells_transformed_filtered_source.dtype.names;
	nx,ny = 2,3
	for i, name in enumerate(names):
	    plt.subplot(nx, ny, i+1)
	    plt.hist(cells_transformed_filtered_source[name]);
	    plt.title(name)
	    if name == 'intensity':
	        plt.xlabel('counts')
	    elif name == 'region':
	        plt.xlabel('region_id')
	    else:
	        plt.xlabel('voxels')
	plt.tight_layout()

	savename=os.path.join(savedir,"transformed_cells_distributions_sizefiltered.png")
	plt.savefig(savename,format='png')
	print(f"Saved {savename}")
	sys.stdout.flush()

	""" Next plots: cells overlaid on atlas """
	pma_file = '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels_60um_edge_80um_vent_erosion.tif'
	pma_vol = tifffile.imread(pma_file)
	xshape,yshape,zshape = pma_vol.shape
	z_planes = list(range(0,zshape,50))
	# z_planes=[200]
	for zplane in z_planes:
		zstr = f'Z{str(zplane).zfill(4)}'

		coordinates = np.hstack([cells_transformed_filtered_source[c][:,None] for c in 'xyz']);
		zplane_depth = 3
		minplane = max(zplane-zplane_depth,0)
		maxplane = zplane+zplane_depth
		zplane_range = range(minplane,maxplane)
		this_plane_coords = np.array([coord for coord in coordinates if coord[-1] in zplane_range])
		try:
			xs = this_plane_coords[:,0]
			ys = this_plane_coords[:,1]
		except:
			xs = []
			ys = []
		# For each quadrant of this slice make a figure
		fig,axes = plt.subplots(figsize=(8,5),nrows=1,ncols=2,sharex=True,sharey=True)
		ax_tissue = axes[0]
		im = pma_vol[zplane]
		mean = im.mean()
		std = im.std()
		ax_tissue.imshow(im,vmin=0,vmax=mean+std*3,cmap='viridis')
		ax_both=axes[1]
		ax_both.imshow(im,vmin=0,vmax=mean+std*3,cmap='viridis')
		ax_both.scatter(xs,ys,s=50,facecolors='none',edgecolors='r')
		plt.suptitle(f"Registered, filtered cells {zstr} shown on atlas, Z{str(zplane_range)}",)
		savename = os.path.join(savedir,f"registered_filtered_cells_overlay_z{zplane}.png")
		plt.savefig(savename,format='png')
		plt.clf() # free up memory
		print(f"Saved {savename}")
		sys.stdout.flush()

	# Finally, make the histograms of counts and count densities in brain regions
	cell_count_csv_file = os.path.join(project_dir,"region_cell_counts_filtered_20px.csv")
	df = pd.read_csv(cell_count_csv_file)
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	h=ax.hist(df.iloc[0,:],bins=50)
	ax.set_yscale('log')
	ax.set_ylabel('Number of regions')
	ax.set_xlabel("Counts")
	title=ax.set_title("Brain region counts histogram (eroded and size-filtered)")
	savename = os.path.join(savedir,"filtered_counts_histogram.png")
	plt.savefig(savename,format='png')
	print(f"Saved {savename}")

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	h=ax.hist(df.iloc[1,:],bins=50)
	ax.set_yscale('log')
	ax.set_ylabel('Number of regions')
	ax.set_xlabel("Count density (per mm^3)")
	title=ax.set_title("Brain region count density histogram (eroded and size-filtered)")
	savename = os.path.join(savedir,"filtered_count_densities_histogram.png")
	plt.savefig(savename,format='png')
	print(f"Saved {savename}")


