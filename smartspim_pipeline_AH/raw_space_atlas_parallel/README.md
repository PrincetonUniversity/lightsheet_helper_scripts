# raw_space_atlas_parallel
This directory contains the pipeline for creating a raw space atlas, using parallel processing (only for step 2).

## How to make the raw space annotation atlas for SmartSPIM 3.6x images

The raw space atlas is useful for cases when you want to look at your data in its full resolution but with the atlas regions overlaid. A common use case for this is c-fos, where it is important to be able to see cells at the highest resolution possible with the atlas boundaries. In normal registration, we register data to the atlas and it gets downsized significantly because the atlases we use are lower resolution (typically 20-25 micron resolution isotropic) than the raw data (which range from 1-5 micron resolution). The trick here is that we are registering the atlas to the data, i.e. going in reverse.

Prerequisites for this are having run the elastix inverse transformations, i.e. the two transformations:
1. resampled regch 488 (moving) -> resampled cellch 642 (fixed)
2. atlas (moving) -> resampled regch 488 (moving)

Where the resampling (aka downsizing) of both regch and cellch is assumed to have been done ahead of time.

## Steps to make the raw space atlas

**Files you need to modify:**

`run_transform_pipeline.sh`

At the top of this file you need to set the following variables:
```
raw_dir=/jukebox/wang/Jess/lightsheet_output/202010_cfos/processed/an001/full_sizedatafld/072420_jv_ymazelearn_an1_1_3x_647_008na_1hfds_z10um_50msec_ch00
# Define the path to where elastix inverse transformations were run
elastix_atlas_to_auto_dir=/jukebox/wang/Jess/lightsheet_output/202010_cfos/processed/an001/ClearMapClusterOutput/elastix_auto_to_atlas
elastix_auto_to_cell_dir=/jukebox/wang/Jess/lightsheet_output/202010_cfos/processed/an001/ClearMapClusterOutput/elastix_cfos_to_auto
# Define where you want your aligned raw atlas files to live
output_dir=/jukebox/wang/ahoag/test_raw_atlas
```

**Files that are optional to modify:**

`transform_annotations_to_fullsize_cfos.py`

By default, the script will use the Princeton Mouse atlas with the hierarchical labeling (i.e. parents always have larger segment IDs than their children). If you want to use a different atlas you will need to modify the line:

```
ann = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels.tif" # Princeton Mouse Annotation Atlas, 16 bit
```

You should now be ready to run the code:

`./run_transform_pipeline.sh`

This bash script submits two sbatch scripts: `slurm_scripts/transform_step1.sh` and `slurm_scripts/transform_step2.sh`, which each call the python script using different command line arguments: `transform_annotations_to_fullsize_cfos.py`

## What actually happens in this script

Basically this script just reverses the steps of registration. The *extremely* confusing part is that in actuality it reverses the steps of inverse registration. To make the distinction between regular registration and inverse registration consider the regular registration of a raw cell channel volume to atlas space. In that scenario we perform these steps:

1. Downsize (aka resample) raw cell channel and raw registration channel data to a shape that is slightly larger (typically 1.3 or 1.4 times) than the atlas in x,y, and z dimensions. Both channels should have the same shape, e.g. (1.3*x_atlas,1.3*y_atlas,1.3*z_atlas) after this step. This both speeds up the registration and can make it more accurate.
2. Reorient the downsized data for both channels so that they are in the sagittal orientation (i.e where z slices are sagittal cuts). Typically the raw data are in the horizontal orientation (i.e. z slices are horizontal cuts), so this amounts to just swapping the x and z axes of the resampled image.
3. Run elastix with the downsized, reoriented cell channel as the moving file and the downsized, reoriented registration channel as the fixed channel. The result is a downsized, reoriented cell channel that has been aligned to the downsized, reoriented cell channel registration channel.
4. Run elastix with the downsized, reoriented registration channel as the moving channel and the atlas tissue volume as the fixed channel. The result is a registration channel that has been aligned to the atlas.

The steps above are for the normal registration. For inverse registration we want to swap the moving and fixed channels in step 3 and 4, like so:
1. (same as step 1 in normal registration)
2. (same as step 2 in normal registration)  
3. Run elastix with the downsized, reoriented cell channel as the **fixed** file and the downsized, reoriented registration channel as the **moving** channel. The result is a downsized, reoriented registration channel that has been aligned to the downsized, reoriented cell channel.
4. Run elastix with the downsized, reoriented registration channel as the **fixed** channel and the atlas tissue volume as the **moving** channel. The result is the atlas which has been aligned to the downsized, reoriented registration channel.

In steps 3 and 4, two files are generated in each step called:
```
TransformParameters.0.txt
TransformParameters.1.txt
```
These files are matrices that represent the transformations between the moving and fixed channels. There are two files because there are actually two different types of transformations that are performed for each step. The 0 file refers to an affine transformation, i.e. a linear transformation with translation. The 1 file refers to a bspline transformation, i.e. a nonlinear transformation. The files `result.0.tif` and `result.1.tif` in the elastix directories correspond to the transformed image (i.e. the "moved" image) after the affine and bspline transformations, respectively. Note that the affine and bspline transformations are entirely independent and typically we don't actually need to run both since we almost always use bspline because it tends to do much better.

Although the python script `transform_annotations_to_fullsize_cfos.py` looks like it uses all 4 of the `TransformParameters.*.txt` files, in actuality we only need the bspline transformation from atlas -> reg, i.e. the `TransformParameters.1.txt` file from the inverse transformation between atlas and reg. You can see this in the line where we call transformix, the program that takes the output of elastix, i.e. the `TransformParameters.*.txt` and uses them to do align images together:  
```
transformix_command_line_call(ann, aldst, transformfiles[-1])
```
The `transformfiles[-1]` there indicates that we're only using the last file in the list `transformfiles`, which is the file I mentioned. This line is the first step in making the raw atlas. What it does is take the atlas annotation volume `ann` and warp it using the bspline inverse transformation that we found between the atlas and the registration channel. The result is the atlas in the space of the sagittal-orientation, downsized registration channel.

So two more things need to happen to get the atlas into the original raw space. First, the sagittal orientation needs to be changed to horizontal orientation. The second is that the downsizing needs to be reversed (aka "zoomed out") so that the atlas has the same dimensions as the true raw data.

The reorientation happens in this line:
```
bigdvann = np.swapaxes(zoom(tann, [1,1,dv0/float(dv1)], order=0),0,2)
```
And also the dorsal-ventral dimension (x dimension in the sagittally oriented volume, z dimension of the resulting reoriented volume) is zoomed out to the raw space dimensions. The resulting volume `bigdvann` is now horizontally oriented and only its x and y dimensions need to be zoomed out. That is done in the `process_slice(z)` function in parallel.


## Something to note
The raw atlas that this script creates is actually aligned to the raw registration channel. However, we often want to show the raw cell channel and the cell annotations with the atlas overlaid. Technically there should be an additional transformation between reg ch -> cell ch after the intitial transformation between atlas -> reg ch. However, the reg channel and cell channel are often essentially already aligned due to how we image the same volume at multiple wavelengths simultaneously. As a result, you won't notice much of a difference if you add in this extra transformation in almost all cases. If it is NOT the case that cell channel and registration channel were imaged at the same resolution together, or you suspect differences between the dimensions of these volumes, then this additional step needs to be run. All that is needed in that case is an extra call to transformix using the correct `TransformParameters.*.txt` file after the initial transformix call. 
