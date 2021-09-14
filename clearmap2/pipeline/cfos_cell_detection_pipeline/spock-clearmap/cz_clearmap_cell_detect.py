# General imports
import os,sys,pickle,time,subprocess,random,shutil
import numpy as np

# ClearMap2 imports
sys.path.append('/jukebox/witten/Chris/python/ClearMap2-master')
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ParallelProcessing.BlockProcessing as bp
import ClearMap.Utils.HierarchicalDict as hdict
import ClearMap.ImageProcessing.Experts.Cells as cells

def process_block(block,savedir,params,verbose):
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
    block_index = block.index[-1]
    block_result = cells.detect_cells_block(block, parameter=params,verbose=verbose)
    block_savename = os.path.join(savedir,f'cells_block{block_index}.p')
    with open(block_savename,'wb') as pkl:
        pickle.dump(block_result,pkl)
    print(f"Saved {block_savename}")
    return block_result

if __name__ == '__main__':
    output_rootpath = '/jukebox/wang/ahoag/for_cz/clearmap2_test_output'
    sample_dir = sys.argv[1].strip().rstrip("/")
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    request_name,sample_name = sample_dir.split('/')[-2:]
    dst_dir = os.path.join(output_rootpath,request_name,sample_name,
        'imaging_request_1/rawdata/resolution_3.6x')
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
    hdict.pprint(cell_detection_parameter)
    print()

    # Split into blocks
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

    block = blocks[jobid]
    print(f"Running cell detection on block {jobid}")
    sys.stdout.flush()
    block_result = process_block(block,savedir=result_dir,
        params=cell_detection_parameter,verbose=True)
    print("Done running cell detection on this block")
    sys.stdout.flush()


#     block_result_list = []
#     print()
#     print('Merging block results into a single data file...')
#     for block in blocks:
#         block_index = block.index[-1]
#         #print(f"Working on block {block_index}")
#         block_savename = os.path.join(directory,f'cells_blocks',f'cells_block{block_index}.p')
#         with open(block_savename,'rb') as pkl:
#             block_result = pickle.load(pkl)
#             block_result_list.append(block_result)
#             header = ['x','y','z']
#             dtypes = [int, int, int]
#             if cell_detection_parameter['shape_detection'] is not None:
#                 header += ['size']
#                 dtypes += [int]
#             measures = cell_detection_parameter['intensity_detection']['measure']
#             header +=  measures
#             dtypes += [float] * len(measures)
#     final_results = np.vstack([np.hstack(r) for r in block_result_list])
#     dt = {'names' : header, 'formats' : dtypes}
#     cells_allblocks = np.zeros(len(final_results), dtype=dt)
#     for i,h in enumerate(header):
#         cells_allblocks[h] = final_results[:,i]
#     savename = ws.filename('cells',postfix='raw')
#     io.write(savename,cells_allblocks)
# print('Saving raw cell detection results to:\n/jukebox/witten/Chris/data/clearmap2/' + fpath + '/cells_raw.npy')
# print()
