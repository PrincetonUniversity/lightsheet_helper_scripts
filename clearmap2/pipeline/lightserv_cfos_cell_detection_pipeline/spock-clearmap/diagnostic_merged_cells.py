# General imports
import os,sys
import matplotlib.pyplot as plt
import numpy as np
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
	savedir = os.path.join(diagnostic_plot_dir,"merged_cells")
	os.makedirs(savedir,exist_ok=True)

	ws = wsp.Workspace('CellMap',directory=project_dir)
	ws.debug = False
	im_source = ws.source('stitched')

	cells_source = ws.source('cells', postfix='raw')

	""" First plot: the histograms of coordinates, size, intensity """
	plt.figure(1); plt.clf();
	names = cells_source.dtype.names;
	nx,ny = 2,3
	for i, name in enumerate(names):
		plt.subplot(nx, ny, i+1)
		plt.hist(cells_source[name]);
		plt.title(name)
		if name in ['source','background']:
			plt.xlabel('counts')
		else:
			plt.xlabel('voxels')

	plt.tight_layout()
	savename=os.path.join(savedir,"raw_cells_distributions.png")
	plt.savefig(savename,format='png')
	print(f"Saved {savename}")
	sys.stdout.flush()
	""" Next plots: cells overlaid on tissue """
	xshape,yshape,zshape = im_source.shape
	z_planes = list(range(0,zshape,250))
	# z_planes=[750]
	for zplane in z_planes:
		zstr = f'Z{str(zplane).zfill(4)}'

		coordinates = np.hstack([cells_source[c][:,None] for c in 'xyz']);
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
		halfy = int(round(xshape/2.)) # intentionally swapped x and y for display purposes below
		halfx = int(round(yshape/2.)) # intentionally swapped x and y for display purposes below
		for xmin in [0,halfx]:
			for ymin in [0,halfy]:
				xmax=xmin+halfx
				ymax=ymin+halfy
				fig,axes = plt.subplots(figsize=(12,5),nrows=1,ncols=2,sharex=True,sharey=True)
				ax_tissue = axes[0]
				im = im_source[:,:,zplane]
				mean = im.mean()
				std = im.std()
				ax_tissue.imshow(im,vmin=0,vmax=mean+std*3,cmap='viridis')
				ax_both=axes[1]
				ax_both.imshow(im,vmin=0,vmax=mean+std*3,cmap='viridis')
				ax_both.scatter(ys,xs,s=50,facecolors='none',edgecolors='r')
				ax_both.set_xlim([xmin,xmax])
				ax_both.set_ylim([ymin,ymax])
				plt.suptitle(f"Unfiltered cells at z={zplane} projected over: Z{str(zplane_range)}, X: ({xmin},{xmax}), Y: ({ymin},{ymax})",)
				savename = os.path.join(savedir,f"raw_cells_overlay_z{zplane}_x{xmin}-{xmax}_y{ymin}-{ymax}.png")
				plt.savefig(savename,format='png')
				plt.clf() # free up memory
				print(f"Saved {savename}")
				sys.stdout.flush()


				


