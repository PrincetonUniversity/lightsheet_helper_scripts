#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 15                 # time (minutes)
#SBATCH -o logs/debug_cell_detect_%A_%a.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/debug_cell_detect_%A_%a.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100gbs 

# tmpdir=/tmp/${UID}_clearmap2_block${SLURM_ARRAY_TASK_ID}
# echo "tmp dir is: $tmpdir"
echo ${SLURM_JOB_ID}
errordir=./errorfiles
authdir=./authfiles
errorfile=${errordir}/errorfile_${SLURM_JOB_ID}
authfile=${authdir}/authfile_${SLURM_JOB_ID}
export XDG_RUNTIME_DIR=$authdir

module load anacondapy/2020.11
conda activate ClearMap
# displayport=$RANDOM
# echo "Using display server: $displayport"
# Sleep for 2*array_job_id seconds to avoid collisions between array jobs
tsleep=`echo "$SLURM_ARRAY_TASK_ID*2+1" | bc`
echo "Sleeping for $tsleep seconds"
sleep $tsleep
# xvfb-run -n $displayport -d -e $errorfile -f $authfile python spock-clearmap/cz_clearmap_cell_detect_debug.py ${sample_dir} ${blocks_per_job} ${output_rootpath}
python spock-clearmap/clearmap_cell_detect_debug.py ${sample_dir} ${blocks_per_job} ${output_rootpath}
EXIT=$?
echo "Exit status is $EXIT"
exit 0
# echo "Removing tmp dir"
# rm -rf $tmpdir
# if [[ "${EXIT}" -eq 135 ]]; then
# 	echo "Exempting exit status 135; instead exiting with 0"
# 	exit 0
# else
# 	exit $EXIT
# fi