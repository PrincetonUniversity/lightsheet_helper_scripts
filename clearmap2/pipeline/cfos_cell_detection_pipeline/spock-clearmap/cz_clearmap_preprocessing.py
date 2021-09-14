# General imports
import os,sys,glob,shutil,json

# ClearMap2 imports
sys.path.append('/jukebox/braininit/lightsheet/ClearMap2')
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ParallelProcessing.BlockProcessing as bp

def fast_scandir(dirname):
	""" gets all folders recursively """
	subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
	for dirname in list(subfolders):
		subfolders.extend(fast_scandir(dirname))
	return subfolders

# Select datasets to analyze
if __name__ == "__main__":
	sample_dir = sys.argv[1].strip().rstrip("/")
	output_rootpath = '/jukebox/wang/ahoag/for_cz/clearmap2_test_output'

	request_name,sample_name = sample_dir.split('/')[-2:]
	src_dir = os.path.join(sample_dir,
		'imaging_request_1/rawdata/resolution_3.6x')
	dst_dir = os.path.join(output_rootpath,request_name,
		sample_name,'imaging_request_1/rawdata/resolution_3.6x')
	if not os.path.exists(dst_dir):
		os.makedirs(dst_dir)
		print(f"Creating dst dir: {dst_dir}")

	# Link and rename image files

	try:
		src_642 = fast_scandir(os.path.join(src_dir,'Ex_642_Em_2_corrected'))[-1]
	except:
		src_642 = os.path.join(src_dir,'Ex_642_Em_2_corrected')

	dst_642 = os.path.join(dst_dir,'Ex_642_Em_2_corrected')

	if not os.path.exists(dst_642):
		print(f"Creating 642 dst dir: {dst_642}")
		os.mkdir(dst_642)
		
	src_files_642 = sorted(glob.glob(src_642 + '/*tif'))
	print()
	print("Sym linking Ch 642 files if not done already")
	for ii,src in enumerate(src_files_642):
		dst_basename = 'Z'+f'{ii}'.zfill(4)+'.tif'
		dst = os.path.join(dst_642,dst_basename)
		if not os.path.exists(dst):
			os.symlink(src,dst)
	try:
		src_488 = fast_scandir(os.path.join(src_dir,'Ex_488_Em_0_corrected'))[-1]
	except:
		src_488 = os.path.join(src_dir,'Ex_488_Em_0_corrected')

	dst_488 = os.path.join(dst_dir,'Ex_488_Em_0_corrected')
	if not os.path.exists(dst_488):
		print(f"Creating 488 dst dir: {dst_488}")
		os.mkdir(dst_488)

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