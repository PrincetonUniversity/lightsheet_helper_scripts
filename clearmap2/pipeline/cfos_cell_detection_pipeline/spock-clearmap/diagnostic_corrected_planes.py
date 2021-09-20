# General imports
import os,sys,glob
import tifffile
import matplotlib.pyplot as plt
from textwrap import wrap

if __name__ == "__main__":
	n_cores = os.cpu_count()
	sample_dir = sys.argv[1].strip().rstrip("/")

	output_rootpath = sys.argv[2].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	project_dir =  os.path.join(output_rootpath,
		request_name,sample_name,"imaging_request_1/rawdata/resolution_3.6x")
	
	ch488_corrected_dir = os.path.join(project_dir,"Ex_488_Em_0_corrected")
	ch488_corrected_files = sorted(glob.glob(ch488_corrected_dir+'/*tif'))
	ch642_corrected_dir = os.path.join(project_dir,"Ex_642_Em_2_corrected")
	ch642_corrected_files = sorted(glob.glob(ch642_corrected_dir+'/*tif'))
	
	channels=["488","642"]
	diagnostic_plot_dir = os.path.join(project_dir,"diagnostic_plots")
	savedir = os.path.join(diagnostic_plot_dir,"corrected_planes")
	os.makedirs(savedir,exist_ok=True)

	for channel in channels:
		corrected_files = eval(f"ch{channel}_corrected_files")
		n_z_planes = len(corrected_files)
		z_planes = list(range(0,n_z_planes,250))
		for z in z_planes:
			z_plane_file = corrected_files[z]
			zstr = os.path.basename(z_plane_file)
			im = tifffile.imread(z_plane_file)
			mean = im.mean()
			std = im.std()
			fig = plt.figure(figsize=(8,5))
			ax=fig.add_subplot(1,1,1)
			ax.imshow(im,vmin=0,vmax=mean+std*3)
			title = f'File: {z_plane_file}'
			ax.set_title("\n".join(wrap(title)),fontsize=10)
			savename = os.path.join(savedir,
				f"channel_{channel}_{zstr}_corrected.png")
			plt.savefig(savename,format='png')
			print(f"Saved {savename}")