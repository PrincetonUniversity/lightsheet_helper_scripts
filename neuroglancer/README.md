# Contents of this folder:
## Notebooks
- brodylab_MRI_atlas_customizations.ipynb
  - A notebook showing how to use Neuroglancer to view an atlas annotation volume with customized labels for the regions.

- detected_cells_cloudvolume.ipynb
  - A notebook showing how to take an array of cell center coordinates (from cell detection output, for example) and convert that so it can be viewable in Neuroglancer.

- detected_cells_cloudvolume (Gaussian filter).ipynb
  - An old notebook in which pixel dilation for detected cells was attempted using a Gaussian filter.

- detected_cells_cloudvolume (Single pixel).ipynb
  - An old notebook in which no pixel dilation was used for detected cells.

- Jess-201904_ymaze_cells.ipynb
  - A notebook in which the clearmap registered cell centers from one of Jess Verpeut's Ymaze experiments were converted using no dilation so that they could be viewable in Neuroglancer (alongside the registered data and atlas)  

- Jess-201904_ymaze_cells_dilate.ipynb
  - A notebook in which the clearmap raw cell centers from one of Jess Verpeut's Ymaze experiments were converted using dilation (so that they were more visible) so that they could be viewable in Neuroglancer (alongside the raw data and raw-space atlas)  

- make_tracing_cells_ng_figure_4NSF.ipynb
  - A notebook in which I (Austin) experimented with various dilation schemes for showing cells from Tom Pisano's tracing data in Neuroglancer, specifically to make a figure for the 2019 NSF equipment proposal (Petaverse).

- neuroglancer_localhost.ipynb
  - A notebook in which I (Austin) experimented with using a locally hosted neuroglancer client (as opposed to Seung lab or Google client) for viewing data.

- preload_segments.ipynb
  - A notebook showing how to preload segments in Neuroglancer when you have a segmentation layer. This is useful for showing specific brain regions or all brain regions. The default is no segments are displayed.

## Scripts
- make_gifs_videos.py
  - Zahra's code for making gifs from Neuroglancer screenshots

- make_precomputed_allenatlas_2017.py
  - How to convert the 2017 Allen brain atlas annotation volume into a precomputed volume that Neuroglancer can view

- make_precomputed_tracing.py
  - Code for converting viral tracing dataset into precomputed format that Neuroglancer can view

- mk_ng_vol.sh
  - Slurm script to run the make_precomputed_tracing.py script on spock/other slurm-based cluster

# Tip: To run most of these notebooks/scripts make an anaconda environment:

conda create --name ng python=3.6 <br>
conda activate ng # go into the environment (or source activate if you're using anaconda2)

\# Get numpy

`pip install numpy`

\# install cloudvolume<br>
\# Clone the repo: https://github.com/seung-lab/cloud-volume, then go into the repo directory wherever you put your git stuff, e.g.:

`cd ~/Git/cloud-volume/`

`pip install -e .[all_viewers]`

\# install neuroglancer

`pip install neuroglancer`

\# install igneous - needed in order to mesh (show objects in 3D in neuroglancer) and to downsample (needed for large volumes)

\# Clone the repo: https://github.com/seung-lab/igneous, then go into the repo directory, e.g.:

`cd ~/Git/igneous/`

`pip install -r requirements.txt`

`python setup.py develop` # Make sure to be inside your ng conda environment for this
