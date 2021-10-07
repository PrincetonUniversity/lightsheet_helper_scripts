### cfos_cell_detection_pipeline

Main script:
```
clearmap_pipeline.sh $request_name $imaging_request $sample_name(optional)
```
This is an automated pipeline that performs whole-brain cell detection, registration of  detected cells to a specified atlas reference space, and creation of a pandas-compatible CSV file containing the breakdown of counts in each brain region in the atlas. The pipeline assumes you have your data organized on `/jukebox/LightSheetData/lightserv/` as a normal Lightserv request.  

The cell detector is: [ClearMap2](https://github.com/ChristophKirst/ClearMap2). Cell detection parameters are currently hardcoded in the main script: `clearmap_pipeline.sh` under the variable name `clearmap_params_file`. 

## Before you run the pipeline
- In the same directory as the main script, do:
```bash
$ mkdir authfiles
$ mkdir logs
```
- Make sure that for each sample in your request there is a 

Also:
- Make a conda environment called `ClearMap` with clearmap2 installed.
- Make a conda environment called `lightsheet` with the following packages installed:
```
pip install numpy scipy tifffile opencv-python 
```
- update the lines at the top of each python file in `spock-clearmap/`,`downsizing/`, and `registration/` that start with `sys.path.append` to the location where you have cloned the light-sheet processing pipeline [BrainPipe](https://github.com/BrainCOGS/BrainPipe)
- Update the variables at the top of `clearmap_pipeline.sh`:
    -`username`: your netid
    -`output_rootpath`: the root location where you want to save the results
    -`clearmap_params_file`: the clearmap parameter pickle file containing the cell detection parameters