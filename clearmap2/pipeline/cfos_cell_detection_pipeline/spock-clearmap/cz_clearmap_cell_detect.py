# General imports
import os,sys,pickle

# ClearMap2 imports
sys.path.append('/jukebox/braininit/lightsheet/ClearMap2')
import ClearMap.IO.Workspace as wsp
import ClearMap.ParallelProcessing.BlockProcessing as bp
import ClearMap.Utils.HierarchicalDict as hdict
import ClearMap.ImageProcessing.Experts.Cells as cells

from functools import partial    
from concurrent.futures import ProcessPoolExecutor

def process_block(savedir,params,block_index):
    """
    ---PURPOSE---
    A function that takes a block (chunk of volume) as input and 
    runs the cells.detect_cells_block() function on it.
    
    We then save the results in an array so that we 
    can just load them later when we want to merge this all together
    ---INPUT---
    block                       A processing block created from bp.split_into_blocks()
    savedir                     The directory in which to save the detected cells from this block
    params                      The cell detection parameter dictionary 
                                that you feed into detect_cells_block()
    verbose                     True or False
    ---OUTPUT---
    block_result      The tuple containing the cell coordinates, shape, intensities 
    It also saves this block_result as a file in your savedir called:
                      "cells_block{block_index}.p" where block_index is ranges from 0 to the number of blocks-1
    """
    # block_index = block.index[-1]
    print("Processing block: ", block_index)
    block = blocks[block_index]
    sys.stdout.flush()
    block_result = cells.detect_cells_block(block, parameter=params)
    block_savename = os.path.join(savedir,f'cells_block{block_index}.p')
    with open(block_savename,'wb') as pkl:
        pickle.dump(block_result,pkl)
    print(f"Saved {block_savename}")
    return "success"

if __name__ == '__main__':
    n_cores = os.cpu_count()

    sample_dir = sys.argv[1].strip().rstrip("/")
    imaging_request = sys.argv[2].strip().rstrip("/")
    output_rootpath = sys.argv[3].strip().rstrip("/")
    blocks_per_job = int(sys.argv[4])

    array_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    request_name,sample_name = sample_dir.split('/')[-2:]
    dst_dir = os.path.join(output_rootpath,request_name,sample_name,
        imaging_request,'rawdata/resolution_3.6x')
    # Initialize ClearMap2 workspace object
   
    ws = wsp.Workspace('CellMap',directory=dst_dir)
    ws.debug = False
    print()
    ws.info()

    # Load ClearMap2 cell detection parameters
    fname = '/jukebox/witten/Chris/data/clearmap2/utilities/cell_detection_parameter.p'
    with open(fname,'rb') as f:
        cell_detection_parameter = pickle.load(f)
    # print('path                    : ' + fname)
    cell_detection_parameter['verbose']=True
    hdict.pprint(cell_detection_parameter)
    print()

    # Get the blocks
    blocks = bp.split_into_blocks(ws.source('stitched'),
                    processes='serial',
                    axes=[2],
                    size_min=10,
                    size_max=30,
                    overlap=5,
                    verbose=True)
    # create output directory if it does not already exist
    result_dir = os.path.join(dst_dir,'cells_blocks')
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    block_mindex = array_id*blocks_per_job
    block_maxdex = min(block_mindex + blocks_per_job,len(blocks))
    block_indices = list(range(block_mindex,block_maxdex))
    # blocks_this_job = blocks[block_mindex:block_maxdex]

    # Process 1 block per core 
    print("Running cell detection on blocks: ", block_indices)
    sys.stdout.flush()
    process_block_par = partial(process_block,
        result_dir,cell_detection_parameter)
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        for block_result in executor.map(process_block_par,block_indices):
            try:
                print(block_result)
                sys.stdout.flush()
            except Exception as exc:
                print(f'generated an exception: {exc}')

    # for block in blocks_thisjob:
    #     print(f"Running cell detection on block {block.index[-1]}")
    #     sys.stdout.flush()
    #     block_result = process_block(block,savedir=result_dir,
    #         params=cell_detection_parameter)
    #     print("Done running cell detection on this block")
    #     sys.stdout.flush()
