#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 2                      # number of cores
#SBATCH -t 30                 # time (minutes)
#SBATCH -o logs/clearmap_cell_detect_%A_%a.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/clearmap_cell_detect_%A_%a.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100gbs 

tmpdir=/tmp/${UID}_clearmap2_block${SLURM_ARRAY_TASK_ID}
echo "tmp dir is: $tmpdir"
if [ ! -d "$tmpdir" ]; then
  echo "Creating tmp dir: $tmpdir"
  mkdir $tmpdir
  chmod 700 $tmpdir
  ls -ld $tmpdir
fi
export XDG_RUNTIME_DIR=$tmpdir

module load anacondapy/2020.11
conda activate ClearMap

# Sleep for time between 0 and 20 seconds depending on array job id. 
# This is to avoid collisions between array jobs
tsleep=`echo "$SLURM_ARRAY_TASK_ID*2" | bc`
echo "Sleeping for $tsleep seconds"
sleep $tsleep
xvfb-run -d python spock-clearmap/cz_clearmap_cell_detect.py ${sample_dir} ${blocks_per_job} ${output_rootpath}
echo "Removing tmp dir"
rm -rf $tmpdir