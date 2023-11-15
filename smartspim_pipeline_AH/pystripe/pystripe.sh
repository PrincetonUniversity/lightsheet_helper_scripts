#!/bin/env bash
#
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=20G
#SBATCH --contiguous
#SBATCH --time=9:00:00
#SBATCH -o logs/spim_pystripe_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/spim_pystripe_%j.err        # STDERR #add _%a to see each array job

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

#required
module load anacondapy/2020.11
conda activate lightsheet

echo "Input directory (path to stitched images):" "$1"
echo "Output directory (path to destriped images):" "$2"

pystripe -i "$1" -o "$2" -s1 256 -s2 512 -w 'db3'
