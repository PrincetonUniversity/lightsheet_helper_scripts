# General imports
import os,sys,pickle,time,subprocess,random,shutil
import numpy as np

# ClearMap2 imports
sys.path.append('/jukebox/witten/Chris/python/ClearMap2-master')
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ParallelProcessing.BlockProcessing as bp
import ClearMap.Utils.HierarchicalDict as hdict

# Select dataset to analyze
fpath = os.path.join(sys.argv[1],'rawdata/resolution_3.6x')

if __name__ == '__main__':

    # Initialize ClearMap2 workspace object
    directory = os.path.join('/jukebox/witten/Chris/data/clearmap2',fpath)
    expression_raw = '/Ex_642_Em_2_corrected/Z<Z,4>.tif'
    expression_auto = '/Ex_488_Em_0_corrected/Z<Z,4>.tif'
    ws = wsp.Workspace('CellMap',directory=directory)
    ws.update(raw=expression_raw)
    ws.update(autofluorescence=expression_auto)
    ws.debug = False
    print()
    ws.info()

    # Load ClearMap2 cell detection parameters
    fname = '/jukebox/witten/Chris/data/clearmap2/utilities/cell_detection_parameter.p'
    with open(fname,'rb') as f:
        cell_detection_parameter = pickle.load(f)
    print('path                    : ' + fname)
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

    # create output directory
    output_dir = os.path.join(directory,'cells_blocks')
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    num_blocks_left = len(blocks)
    blocks_left = np.arange(0, num_blocks_left)
    batch_size = 10
    block_iters = {i : 0 for i in blocks_left}
    iters = 0
    max_iter = 100
    while(len(blocks_left) > 0):
        if (iters > max_iter):
            break
        print()
        print("Iteration " + str(iters) + ": " + str(len(blocks_left)) + " blocks remaining")

        if len(blocks_left) < batch_size:
            batch_size = len(blocks_left)
        blocks_to_run = random.sample(list(blocks_left),batch_size)
        blocks_to_run.sort()

        blocks_to_run_str = str(blocks_to_run)
        blocks_to_run_str = blocks_to_run_str.replace(" ", "")
        blocks_to_run_str = blocks_to_run_str.replace("[", "")
        blocks_to_run_str = blocks_to_run_str.replace("]", "")

        slurmpath = '/jukebox/witten/Chris/python/spock-clearmap/slurm/clearmap_blocks_batch.sh'
        run_job = slurmpath + ' ' + blocks_to_run_str + ' ' + fpath
        jobid = subprocess.check_output(run_job, shell=True)
        jobid = int(jobid)

        time.sleep(60) #in seconds

        finished_blocks = []
        while(len(blocks_to_run) > 0):
            blockid = blocks_to_run[0]
            check_status = 'sacct -X -j ' + str(jobid) + '_' + str(blockid)
            status = subprocess.check_output(check_status, shell=True)
            status = str(status)
            if ('PENDING' in status):
                print("Block pending: " + str(blockid))
                time.sleep(60)
                continue
            elif ('RUNNING' in status):
                print("Block running: " + str(blockid))
                time.sleep(60)
                continue
            elif ('COMPLETED' in status):
                block_savename = os.path.join(directory,'cells_blocks',f'cells_block{blockid}.p')
                if (os.path.isfile(block_savename)):
                    finished_blocks.append(blockid)
                    blocks_to_run.remove(blockid)
                    print("Block saved: " + str(blockid))
                else:
                    print("Block failed to save: " + str(blockid))
                    blocks_to_run.remove(blockid)
            elif ('FAILED' in status):
                block_iters[blockid] = block_iters[blockid] + 1
                blocks_to_run.remove(blockid)
                print("Block failed to run: " + str(blockid))
            else:
                print("Block submitted: " + str(blockid))
                time.sleep(60)
                continue

        blocks_left = np.delete(blocks_left, np.isin(blocks_left,finished_blocks))
        iters = iters + 1

    block_result_list = []
    print()
    print('Merging block results into a single data file...')
    for block in blocks:
        block_index = block.index[-1]
        #print(f"Working on block {block_index}")
        block_savename = os.path.join(directory,f'cells_blocks',f'cells_block{block_index}.p')
        with open(block_savename,'rb') as pkl:
            block_result = pickle.load(pkl)
            block_result_list.append(block_result)
            header = ['x','y','z']
            dtypes = [int, int, int]
            if cell_detection_parameter['shape_detection'] is not None:
                header += ['size']
                dtypes += [int]
            measures = cell_detection_parameter['intensity_detection']['measure']
            header +=  measures
            dtypes += [float] * len(measures)
    final_results = np.vstack([np.hstack(r) for r in block_result_list])
    dt = {'names' : header, 'formats' : dtypes}
    cells_allblocks = np.zeros(len(final_results), dtype=dt)
    for i,h in enumerate(header):
        cells_allblocks[h] = final_results[:,i]
    savename = ws.filename('cells',postfix='raw')
    io.write(savename,cells_allblocks)
print('Saving raw cell detection results to:\n/jukebox/witten/Chris/data/clearmap2/' + fpath + '/cells_raw.npy')
print()
