# General imports
import os,sys,glob,shutil,json
cwd = os.getcwd()
utils_fld = os.path.join(cwd,"utils")
sys.path.append(utils_fld)
from pipeline_utils import fast_scandir
# ClearMap2 imports
sys.path.append('/jukebox/braininit/lightsheet/ClearMap2')
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ParallelProcessing.BlockProcessing as bp

# Select datasets to analyze
if __name__ == "__main__":
	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	src_dir = os.path.join(sample_dir,imaging_request,
		'rawdata/resolution_3.6x')
	print(src_dir)
	dst_dir = os.path.join(output_rootpath,request_name,
		sample_name,imaging_request,'rawdata/resolution_3.6x')
	
	if not os.path.exists(dst_dir):
		os.makedirs(dst_dir)
		print(f"Creating dst dir: {dst_dir}")

	# Link and rename image files. Need to check if files are in old 
	# old folder convention e.g. Ex_642_Em_2
	# or new folder convention e.g. Ex_647_Em_680
	old_src_642 = os.path.join(src_dir,'Ex_642_Em_2_corrected')
	#new_src_642 = os.path.join(src_dir,'Ex_647_Em_690_stitched')
	new_src_642 = os.path.join(src_dir,'Ex_639_Em_690_stitched')
	src_642  = ""
	if os.path.exists(old_src_642):
		try:
			src_642 = fast_scandir(old_src_642)[-1]
			print(src_642)
		except:
			src_642 = old_src_642
			print(src_642)
	elif os.path.exists(new_src_642):
		try:
			src_642 = fast_scandir(new_src_642)[-1]
			print(src_642)
		except:
			src_642 = new_src_642
			print(src_642)

	print(src_642)
	dst_642 = os.path.join(dst_dir,'Ex_642_Em_2_corrected')

	os.makedirs(dst_642,exist_ok=True)
		
	src_files_642 = sorted(glob.glob(src_642 + '/*tif'))
	n_src_files_642 = len(src_files_642)
	print(f"have {n_src_files_642} ch642 corrected planes")
	print()
	print("Sym linking Ch 642 files if not done already")
	for ii,src in enumerate(src_files_642):
		dst_basename = 'Z'+f'{ii}'.zfill(4)+'.tif'
		dst = os.path.join(dst_642,dst_basename)
		if not os.path.exists(dst):
			os.symlink(src,dst)
	## Now 488
	src_488 = ""
	old_src_488 = os.path.join(src_dir,'Ex_488_Em_0_corrected')
	new_src_488 = os.path.join(src_dir,'Ex_488_Em_525_stitched')
	if os.path.exists(old_src_488):
		try:
			src_488 = fast_scandir(old_src_488)[-1]
		except:
			src_488 = old_src_488
	elif os.path.exists(new_src_488):
		try:
			src_488 = fast_scandir(new_src_488)[-1]
		except:
			src_488 = new_src_488

	dst_488 = os.path.join(dst_dir,'Ex_488_Em_0_corrected')
	os.makedirs(dst_488,exist_ok=True)

	print("Sym linking Ch 488 files if not done already")
	src_files_488 = sorted(glob.glob(src_488 + '/*tif'))
	for ii,src in enumerate(src_files_488):
		dst_basename = 'Z'+f'{ii}'.zfill(4)+'.tif'
		dst = os.path.join(dst_488,dst_basename)
		if not os.path.exists(dst):
			os.symlink(src,dst)

	# Initialize ClearMap2 workspace object
	# directory = os.path.join(output_path,fpath)
	expression_raw = '/Ex_642_Em_2_corrected/Z<Z,4>.tif'
	ws = wsp.Workspace('CellMap',directory=dst_dir)
	ws.update(raw=expression_raw)
	ws.debug = False
	print()
	ws.info()
	
	# Finally, figure out how many blocks there so we can know how many array jobs to make
	print("Determining number of blocks")
	sys.stdout.flush()
	blocks = bp.split_into_blocks(ws.source('raw'),
                    processes='serial',
                    axes=[2],
                    size_min=10,
                    size_max=30,
                    overlap=5,
                    verbose=True)
	
	n_blocks = len(blocks)
	block_json_file = os.path.join(dst_dir,'block_processing_info.json')
	block_data = {
		'n_blocks':n_blocks
		}
	with open(block_json_file,'w') as outfile:
		json.dump(block_data,outfile)
	print(f"Wrote block info to file: {block_json_file}")
