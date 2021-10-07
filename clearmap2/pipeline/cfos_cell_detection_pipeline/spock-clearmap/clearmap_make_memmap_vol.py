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
	imaging_request = sys.argv[2].strip().rstrip("/")

	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	src_dir = os.path.join(sample_dir,imaging_request,
		'rawdata/resolution_3.6x')
	dst_dir = os.path.join(output_rootpath,request_name,
		sample_name,imaging_request,'rawdata/resolution_3.6x')
	os.makedirs(dst_dir,exist_ok=True)
	print(f"Creating dst dir: {dst_dir}")

	# Link and rename image files

	try:
		src_642 = fast_scandir(os.path.join(src_dir,'Ex_642_Em_2_corrected'))[-1]
	except:
		src_642 = os.path.join(src_dir,'Ex_642_Em_2_corrected')

	dst_642 = os.path.join(dst_dir,'Ex_642_Em_2_corrected')

	os.makedirs(dst_642,exist_ok=True)
		
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

	# Create combined npy volume -- this is CPU heavy
	source = ws.source('raw')
	sink = ws.filename('stitched')
	sys.stdout.flush()
	if not os.path.exists(sink):
		print("Converting raw files into single stitched.npy file")
		sys.stdout.flush()
		io.convert(source,sink,verbose=True)
		print("Successfully created stitched.npy file")
	else:
		print("Stitched.npy file already exists")
	sys.stdout.flush()