Z. Dhanerawala. April 2020.
zahra.dhanerawala@gmail.com.

This directory is a central location for files and documentation related to the Princeton Mouse Atlas used in the Princeton Neuroscience Institute. For original info refer to `wang/pisano/Python/pytom/jupyter_notebooks/registration/finalscripts_lsatlas` or Pisano et al. (2020).

Contents in this folder include:
1. `allen_atlas`: Atlas files relevant to the Allen Mouse Brain Atlas.
2. `sagittal_atlas_20um_iso.tif`: The 3D template Princeton Mouse Atlas used to register light-sheet whole brain volumes. 20um/pixel resolution in X,Y,Z.
3. `cb_sagittal_atlas_20um_iso.tif`: The 3D template Princeton Mouse Atlas used to register light-sheet whole brain volumes, cropped to include only cerebellum. Useful in experiments where the cerebellum is cut off and cleared separately from the whole brain. 20um/pixel resolution in X,Y,Z.
4. `annotation_sagittal_atlas_20um_iso.tif`: The 3D annotation file corresponding to the Princeton Mouse Atlas. 20um/pixel resolution in X,Y,Z.
5. `cb_annotation_sagittal_atlas_20um_iso.tif`: The 3D annotation file corresponding to the Princeton Mouse Atlas, cropped to include only cerebellum. 20um/pixel resolution in X,Y,Z. Used along with `cb_sagittal_atlas_20um_iso.tif`
6. `annotation_sagittal_atlas_20um_iso_16bit.tif`: The 3D annotation file corresponding to the Princeton Mouse Atlas converted to a 16-bit file by reassigning out-of-range integer values to within the range of the 16-bit image. 20um/pixel resolution in X,Y,Z.
7. `annotation_sagittal_atlas_20um_iso_16bit_60um_edge_80um_vntric_erosion.tif`: DEPRECATED: Custom Princeton Mouse Atlas annotation file where edges and ventricles of the brain volumes are eroded by a certain distance. Used in T. Pisano's tracing project.
8. `ls_id_table_w_voxelcounts.xlsx`: The ontology file with an additional column showing total voxels in each structure (row).
9. `ls_id_table_w_voxelcounts_16bit.xlsx`: The ontology dataframe for the corresponding 16-bit annotation volume. Useful for visualization in full-resolution data.
