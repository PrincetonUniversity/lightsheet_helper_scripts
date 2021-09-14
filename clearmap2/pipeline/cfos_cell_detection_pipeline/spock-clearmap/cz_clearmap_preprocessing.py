# General imports
import os,sys,glob,shutil

rootpath = "/jukebox/LightSheetData/lightserv/cz15"
output_path = '/jukebox/witten/Chris/data/clearmap2'
# ClearMap2 imports
sys.path.append('/jukebox/witten/Chris/python/ClearMap2-master')
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io

# Select dataset to analyze
sample_name = sys.argv[1]
fpath = os.path.join(sample_name,'rawdata/resolution_3.6x')

# Link and rename image files
def fast_scandir(dirname):
    """ gets all folders recursively """
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders

src_dir = os.path.join(rootpath,fpath)
dst_dir = os.path.join(output_path,fpath)

try:
    src_642 = fast_scandir(os.path.join(src_dir,'Ex_642_Em_2_corrected'))[-1]
except:
    src_642 = os.path.join(src_dir,'Ex_642_Em_2_corrected')

dst_642 = os.path.join(dst_dir,'Ex_642_Em_2_corrected')
if os.path.isdir(dst_642):
    shutil.rmtree(dst_642)
os.makedirs(dst_642)
src_files = sorted(glob.glob(src_642 + '/*tif'))
for ii,src in enumerate(src_files):
    dst_basename = 'Z'+f'{ii}'.zfill(4)+'.tif'
    dst = os.path.join(dst_642,dst_basename)
    os.symlink(src,dst)
try:
    src_488 = fast_scandir(os.path.join(src_dir,'Ex_488_Em_0_corrected'))[-1]
except:
    src_488 = os.path.join(src_dir,'Ex_488_Em_0_corrected')

dst_488 = os.path.join(dst_dir,'Ex_488_Em_0_corrected')

if os.path.isdir(dst_488):
    shutil.rmtree(dst_488)

os.makedirs(dst_488)
src_files = sorted(glob.glob(src_488 + '/*tif'))
for ii,src in enumerate(src_files):
    dst_basename = 'Z'+f'{ii}'.zfill(4)+'.tif'
    dst = os.path.join(dst_488,dst_basename)
    if not os.path.exists(dst):
        os.symlink(src,dst)

# Initialize ClearMap2 workspace object
directory = os.path.join(output_path,fpath)
expression_raw = '/Ex_642_Em_2_corrected/Z<Z,4>.tif'
ws = wsp.Workspace('CellMap',directory=directory)
ws.update(raw=expression_raw)
ws.debug = False
print()
ws.info()

# Create stitcehd image volume
source = ws.source('raw')
sink = ws.filename('stitched')
io.convert(source,sink,verbose=True)
