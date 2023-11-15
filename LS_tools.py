"""
TO DO:
    -__check_sample_exists function
    -figuring out dorsal/ventral
    -getting rawdata and output paths for dorsal/ventral
    -return list of lasers    
    -implement object refresh to update paths after jobs have been run
    -incorporate params cells from jupyter notebook into this somehow... 
    -just a class with a bunch of constants for each set of user/param pairs
"""

import os

RAWDATA = "rawdata"
OUTPUT = "output"
VIZ = "viz"

SS_RES  = ["resolution_3.6x", "resolution_3.6x_ventral_up", "resolution_15x"]
LV_RES = ["resolution_1.1x", "resolution_1.3x"]

"""
Orientation is always a boolean. True is dorsal and ventral is False.
"""
class LS_tools:
    
    """
    --CONSTRUCTOR--
    Takes in a user, request name, and resolution and finds all relevant samples.
    """
    def __init__(self, user, request_name, resolution, imaging_request=1, 
                 wk_dir="/jukebox/LightSheetData/lightserv", 
                 processing_request=1):
        self.user = user
        self.request_name = request_name
        self.imaging_request = "imaging_request_" + str(imaging_request)
        self.resolution = resolution
        self.ventral_up_resolution = resolution + "_ventral_up"
        self.processing_request = "processing_request_" + str(processing_request)
        self.wk_dir = wk_dir
        self.working_path=os.path.join(self.wk_dir, self.user, self.request_name) 
        if not (os.path.exists(self.working_path)):
            raise ValueError("Working path does not exist")
        
        all_samples = [d for d in os.listdir(self.working_path) if d != ".DS_Store"] 
        samples = [s for s in all_samples if self.imaging_request in os.listdir(os.path.join(self.working_path, s))]
        self.all_samples = sorted(samples)
        if (len(self.all_samples) == 0):
            raise ValueError("No samples found")
        
        self.dorsal_samples = self.__find_samples(self.resolution)
        self.ventral_samples = self.__find_samples(self.ventral_up_resolution)
        
        self.ch_dict_dor = {}
        if resolution in SS_RES:
            for sample in self.dorsal_samples:
                self.ch_dict_dor[sample] = self.__find_sample_channels(sample)
        
        self.ch_dict_ven = {}
        if resolution in SS_RES:
            for sample in self.ventral_samples:
                self.ch_dict_ven[sample] = self.__find_sample_channels(sample)
    
    def get_working_path(self):
        return self.working_path

    def get_all_samples(self):
        return self.all_samples
    
    def get_dorsal_samples(self):
        return self.dorsal_samples
    
    def get_ventral_samples(self):
        return self.ventral_samples

    """
    Add ventral
    """
    def get_sample_channels(self, sample, orientation=True):
        if orientation:
            return self.ch_dict_dor[sample]
        else:
            return self.ch_dict_ven[sample]
    
    def get_sample_channel(self, sample, channel_type, laser, orientation=True):
        sample_dict = self.get_sample_channels(sample, orientation)
        sample_channel = sample_dict[channel_type]
        channel = [ch for ch in sample_channel if str(laser) in ch and not "raw" in ch]

        lasers = [639, 642, 647]

        if (laser in lasers):
            channel = [ch for ch in sample_channel if str(lasers[0]) in ch or 
                                                      str(lasers[1]) in ch or
                                                      str(lasers[2]) in ch]
    

        if len(channel) == 0:
            raise ValueError("Specified channel does not exist for current sample")
        return channel[0]
        

    def get_rawdata_path(self, sample, orientation=True):
        
        #res_type must be ventral or dorsal
 
        if (self.resolution in SS_RES):
           return os.path.join(self.working_path, sample, self.imaging_request, RAWDATA, self.resolution)
        else:
           return os.path.join(self.working_path, sample, self.imaging_request, RAWDATA, self.resolution, sample)

   
    """
    FIX THIS
    """
    def get_output_path(self, sample, res_type):
        self.__check_sample_ind(sample, res_type)
        
        if (self.resolution in LV_RES):
           return os.path.join(self.working_path, sample, self.imaging_request, OUTPUT, self.processing_request,
                               self.resolution) 
    
    def get_viz_path(self, sample, ls=False):
        if (ls):
            return os.path.join(self.working_path, sample, self.imaging_request, VIZ)
        else: 
            return os.path.join("/jukebox/LightSheetData/neuroglancer/public/lightserv", 
                                self.user, self.request_name, sample, self.imaging_request, VIZ)
    
    def get_channel_path(self, sample, channel_type, laser, orientation=True):
        rawdata_sample_path = self.get_rawdata_path(sample)
        channel = self.get_sample_channel(sample, channel_type, laser)
        return os.path.join(rawdata_sample_path, channel)
 
    #working_path: stitched or corrected path
    def get_nested_tif_dir(self, sample, channel_type, laser, orientation=True): 
        output_path = None
        channel_path = self.get_channel_path(sample, channel_type, laser)
        if (len(os.listdir(channel_path)) > 10):
            output_path = channel_path
        else: 
            RES_path = [d.path for d in os.scandir(channel_path) if d.is_dir() and 'RES' in d.path]
            two_up_path = [d.path for d in os.scandir(RES_path[0]) if d.is_dir() and os.path.basename(d.path).isnumeric()]
            one_up_dir = [d.path for d in os.scandir(two_up_path[0]) if d.is_dir()]
            output_path = one_up_dir[0]
        
        return output_path
        
    
    def get_raw_atlas_pipeline_directories(self, sample, channel_type, laser, elastix_dir_name, ra_dir_name):
        raw = self.get_nested_tif_dir(sample, channel_type, laser)
        elastix_atlas_to_auto = os.path.join(self.get_rawdata_path(sample), elastix_dir_name) 
        elastix_auto_to_cell = os.path.join(self.get_rawdata_path(sample), elastix_dir_name, 
                                            "reg_to_cell")
                                            #"Ex_488_Em_525_downsized_to_Ex_647_Em_690_downsized")
                                            #"488_to_647")                                           
        output = os.path.join(self.get_rawdata_path(sample), ra_dir_name)

        return raw, elastix_atlas_to_auto, elastix_auto_to_cell, output
    
    def get_raw_atlas_path(self, sample, ra_dir_name):
        raw_atlas_path = os.path.join(self.get_rawdata_path(sample), ra_dir_name, "transformed_annotations",
                                      "single_tifs")
        
        return raw_atlas_path
    
    def __find_sample_channels(self, sample, orientation=True):
        #channel_type: "raw", "stitched", or "corrected"
        #laser: 0, 1, or 2
        path = os.path.join(self.get_rawdata_path(sample))
        raw_channels = sorted([c for c in os.listdir(path) if not ("stitched" in c) and not ("corrected" in c) 
                               and not ("downsized" in c) and not ("MIP" in c) and "Ex" in c]) 
        stitched_channels = sorted([c for c in os.listdir(path) if "stitched" in c and not "old" in c])
        corrected_channels = sorted([c for c in os.listdir(path) if "corrected" in c and "Ex" in c])
        downsized_channels = sorted([c for c in os.listdir(path) if "downsized" in c]) #add_lasers_check!!
        
        ch_dict = {"raw" : raw_channels, "stitched" : stitched_channels, "corrected" : corrected_channels,
                   "downsized" : downsized_channels}
        return ch_dict
    
    def __find_samples(self, resolution):
        samples = [s for s in self.all_samples if resolution in os.listdir(os.path.join(self.working_path, s, self.imaging_request, RAWDATA))]
        if (len(samples) > 0):
           return samples
        else:
           return []
        
        

        
        


