#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 12                      # number of cores
#SBATCH -t 200                # time (minutes)
#SBATCH -o /scratch/zmd/logs/mk_ngl_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/mk_ngl_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 50000 #50 gbs

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/5.3.1
. activate lightsheet

python make_precomputed_tracing.py 20161207_db_bl6_lob6a_850r_53hr cells
