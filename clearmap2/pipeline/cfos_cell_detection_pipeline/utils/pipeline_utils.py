# pipeline_utils.py
import os

def fast_scandir(dirname):
	""" gets all folders recursively """
	subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
	for dirname in list(subfolders):
		subfolders.extend(fast_scandir(dirname))
	return subfolders