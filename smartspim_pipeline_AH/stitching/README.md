# README

This version of the SmartSPIM stitching pipeline properly catches errors. The only file you need to run is the bash script, `spim_stitching_pipeline.sh`. You pass it the input directory where the files to be stitched live, followed by the output "stitched" directory where you want the stitched files to live like: 
```
./spim_stitching_pipeline.sh /path/to/terastitcher_folder_hierarchy /path/to/stitched
```

This bash script calls the sbatch scripts in the `slurm_scripts` folder which each call the python file `spim_stitch.py`. 

**Important** By default, all of the .out and .err logs from the sbatch scripts will be written to a folder called logs/. This folder must exist in the directory where you end up running  `spim_stitching_pipeline.sh` or else the pipeline will fail at step 0, the import step. Make sure you create this `logs` directory before you run the script. Although the slurm scripts are in a subfolder called `slurm_scripts` the `logs` folder is NOT supposed to be in this folder. If in doubt, use full paths to the logs folder. You can change the output log directory in the `.sh` files in the `slurm_scripts` folder if you prefer. Just make sure that the log directory that you choose exists before you run the pipeline.