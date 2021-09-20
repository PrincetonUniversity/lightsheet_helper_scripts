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
	output_rootpath = sys.argv[2].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	project_dir =  os.path.join(output_rootpath,
		request_name,sample_name,"imaging_request_1/rawdata/resolution_3.6x")

	ws = wsp.Workspace('CellMap',directory=project_dir)
	ws.debug = False
	
	source = ws.source('stitched')

	diagnostic_plot_dir = os.path.join(project_dir,"diagnostic_plots")
	savedir = os.path.join(diagnostic_plot_dir,"memmap_vol")
	os.makedirs(savedir,exist_ok=True)

	source = ws.source('stitched')
	n_z_planes = source.shape[-1]
	z_planes = list(range(0,n_z_planes,250))
	for z in z_planes:
	    zstr = f'Z{str(z).zfill(4)}'
	    im = source[:,:,z]
	    mean = im.mean()
	    std = im.std()
	    fig = plt.figure(figsize=(8,5))
	    ax=fig.add_subplot(1,1,1)
	    ax.imshow(np.fliplr(np.rot90(im,k=3)),vmin=0,vmax=mean+std*3)
	    ax.set_title(f'Memmap volume (stitched.npy) ch 642 z={z}')
	    savename = f"{savedir}/memmap_vol_ch642_{zstr}.png"
	    plt.savefig(savename,format='png')
	    print(f"Saved {savename}")


			


