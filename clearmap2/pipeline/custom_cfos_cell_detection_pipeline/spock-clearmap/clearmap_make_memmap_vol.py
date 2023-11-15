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
	sample_dir = sys.argv[1].strip()
	request_name = sys.argv[2].strip()
	output_rootpath = sys.argv[3].strip().rstrip("/")

	path,sample_name = sample_dir.split('/')[-2:]
	src_dir = sample_dir
	dst_dir = os.path.join(output_rootpath, request_name, sample_name)
	
	#os.makedirs(dst_dir,exist_ok=True)
	#print(f"Creating dst dir: {dst_dir}")

	# Initialize ClearMap2 workspace object
	# directory = os.path.join(output_path,fpath)
	expression_raw = 'Ex_642_Em_2_corrected/Z<Z,4>.tif'
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