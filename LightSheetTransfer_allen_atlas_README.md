Z. Dhanerawala. April 2020.
zahra.dhanerawala@gmail.com.

This directory is a central location for files and documentation related to the Allen Mouse Brain Atlas used in the Princeton Neuroscience Institute.

Contents in this folder include:
1. `allen_documentation`: Documentation regarding atlas and CCF version downloaded from the Allen Brain Institute website.
2. `allen.json`: The original ontology file downloaded from the Allen Brain Institue website.
3. `allen_id_table.xlsx`: The ontology json file converted to a dataframe format for analysis and visualization.
4. `allen_id_table_w_voxel_counts.xlsx`: The ontology dataframe with an additional column showing total voxels in each structure (row). Made by Zahra D. using https://github.com/PrincetonUniversity/lightsheet_helper_scripts/blob/master/registration/find_total_voxel_counts_in_annotations.py.
5. `allen_id_table_w_voxel_counts_totals.xlsx`: The ontology file with an additional column showing total voxels in each structure (row). Made by Kelly S.
6. `allen_id_table_w_voxel_counts_16bit.xlsx`: The ontology dataframe for the corresponding 16-bit annotation volume. Useful for visualization in full-resolution data.
7. `annotation_2017_25um_annotation_template_25_sagittal_forDVscans_75um_erosion.tifsagittal_forDVscans.nrrd`: The 3D annotation file (2017 version) downloaded from the Allen Brain Institue website. 25um/pixel resolution in X,Y,Z.
8. `annotation_2017_25um_sagittal_forDVscans_16bit.tif`: The 3D annotation file (2017 version) converted to a 16-bit file by reassigning out-of-range integer values to within the range of the 16-bit image. 25um/pixel resolution in X,Y,Z.
9. `annotation_2017_25um_sagittal_forDVscans_left_side_only.tif`: The 3D annotation file (2017 version) where the right-side of the brain is zeroed out. Useful to acquire cell counts from the left side of the brain only. 25um/pixel resolution in X,Y,Z.
10. `annotation_2017_25um_sagittal_rank3.tif`: Test file where all structures are summed up to the 3rd level of the Allen Mouse Brain hierarchy. 25um/pixel resolution in X,Y,Z.
11. `annotation_template_25_sagittal_forDVscans.tif`: DEPRECATED. The 3D annotation file (unknown version) downloaded from the Allen Brain Institute website.
12. `annotation_template_25_sagittal_forDVscans_75um_erosion.tif`: DEPRECATED. The 3D annotation file (unknown version) downloaded from the Allen Brain Institue website. The pixels at the edges of the brain volumes were zeroed out by 75 um (3 pixels) using https://github.com/PrincetonUniversity/lightsheet_helper_scripts/blob/master/registration/erode_atlas_tp_version.py
13. `average_template_25_sagittal_forDVscans.tif`: The 3D template Allen Mouse Brain used to register light-sheet whole brain volumes, downloaded from the Allen Brain Institute website. 25um/pixel resolution in X,Y,Z.
14. `average_template_25_sagittal_half.tif`: DEPRECATED: Custom Allen Mouse Brain template used to register to half brain volumes. The dorsal and ventral side is flipped in this sagittal volume as well.
15. `cerebellum_atlas`: 3D template and annotation volumetric images of the cerebellum for cerebellum-to-cerebellum partial registration. Useful in experiments where the cerebellum is cut off and cleared separately from the whole brain.
