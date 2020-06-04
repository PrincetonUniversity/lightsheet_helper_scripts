#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -n 1                      # number of cores
#SBATCH -t 20                 # time (minutes)
#SBATCH -o logs/maincells.out        # STDOUT
#SBATCH -e logs/maincells.err        # STDERR

viz_dir=/jukebox/scratch/zmd/
animal_id=20170116_tp_bl6_lob6b_lpv_07

## Now raw cells
echo "Jobids for raw cells precomputed steps 0, 1 and 2:"
#Step 0: Make info file and layer directory

OUT3=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_cells_step0.sh) 
echo $OUT3

#Step 1: Upload raw cells to vol (writes precomputed data to disk)
OUT4=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_cells_step1.sh) 
echo $OUT4

# Step 2: Make downsamples (higher mips) 
OUT5=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_cells_step2.sh) 
echo $OUT5