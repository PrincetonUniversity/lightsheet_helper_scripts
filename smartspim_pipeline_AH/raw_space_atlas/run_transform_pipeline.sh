#!/bin/env bash
brain=/jukebox/LightSheetData/lightserv/cz15/zimmerman_01/zimmerman_01-001/imaging_request_1/
OUT0=$(sbatch --parsable --export=ALL,brain=${brain} slurm_scripts/transform.sh )
echo $OUT0
