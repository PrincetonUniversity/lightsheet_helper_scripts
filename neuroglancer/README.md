# Install notes for cloudvolume/neuroglancer

conda create --name ng python=3.6
conda activate ng # go into the environment (or source activate if you're using anaconda2)
# Get numpy
pip install numpy
# install cloudvolume
# Clone the repo, then go into the repo directory wherever you put your git stuff, e.g.:
cd ~/Git/cloud-volume/
pip install -e .[all_viewers] 
# install neuroglancer
pip install neuroglancer
# get igneous - needed in order to show objects in 3D in neuroglancer
cd ~/Git/igneous/
pip install -r requirements.txt
python setup.py develop # Make sure to be inside your ng conda environment for this