import os

WK_DIR = "/jukebox/LightSheetData/lightserv/"
RAWDATA = "rawdata"
OUTPUT = "output"
VIZ = "viz"

SS_RES  = ["resolution_3.6x"]
LV_RES = ["resolution_1.1x", "resolution_1.3x"]
    
RAW_CHANNELS = ["Ex_488_Em_0", "Ex_561_Em_1", "Ex_642_Em_2"]
STITCHED_CHANNELS = ["Ex_488_Em_0_stitched", "Ex_561_Em_1_stitched", "Ex_642_Em_2_stitched"]
CORRECTED_CHANNELS = ["Ex_488_Em_0_corrected", "Ex_561_Em_1_corrected", "Ex_642_Em_2_corrected"]
DOWNSIZED_CHANNELS = ["Ex_488_Em_0_downsized", "Ex_561_Em_1_downsized", "Ex_642_Em_2_downsized"]

class LS_tools:
    
    def __init__(self, user, request_name, imaging_request, resolution, processing_request="processing_request_1"):
        self.user = user
        self.request_name = request_name
        self.imaging_request = "imaging_request_" + str(imaging_request)
        self.resolution = resolution
        self.processing_request = processing_request
        self.working_path=os.path.join(WK_DIR, self.user, self.request_name) 
        
        all_samples = [d for d in os.listdir(self.working_path)] 
        samples = [s for s in all_samples if self.imaging_request in os.listdir(os.path.join(self.working_path, s))]
        self.samples = sorted(samples)
    
    def get_working_path(self):
        return self.working_path

    def get_samples(self):
        return self.samples 

    def get_rawdata_sample_path(self, sample_ind, channel_type=None, laser=None):
        self.__check_sample_ind(sample_ind)
 
        sample = self.get_samples()[sample_ind]
        if (self.resolution in SS_RES):
           channel = self.__get_channel_type(channel_type, laser)
           return os.path.join(self.working_path, sample, self.imaging_request, RAWDATA, self.resolution, channel)
        else:
           return os.path.join(self.working_path, sample, self.imaging_request, RAWDATA, self.resolution, sample)
   
    def get_output_sample_path(self, sample_ind):
        self.__check_sample_ind(sample_ind)
        
        if (self.resolution in SS_LV):
           return os.path.join(self.working_path, sample, self.imaging_request, OUTPUT, self.processing_request,
                               self.resolution)  
 
    #working_path: stitched or corrected path
    def get_nested_tif_dir(self, sample_ind, channel_type, laser): 
        output_path = None
        rawdata_sample_path = self.get_rawdata_sample_path(sample_ind, channel_type, laser)
        if (len(os.listdir(rawdata_sample_path)) > 10):
            output_path = rawdata_sample_path
        else: 
            RES_path = [d.path for d in os.scandir(rawdata_sample_path) if d.is_dir() and 'RES' in d.path]
            two_up_path = [d.path for d in os.scandir(RES_path[0]) if d.is_dir() and os.path.basename(d.path).isnumeric()]
            one_up_dir = [d.path for d in os.scandir(two_up_path[0]) if d.is_dir()]
            output_path = one_up_dir[0]
        
        return output_path
    
    def __get_channel_type(self, channel_type, laser):
        #channel_type: "raw", "stitched", or "corrected"
        #laser: 0, 1, or 2
        if (not laser in [0, 1, 2]):
           raise ValueError("Laser must be 0, 1, or 2")
 
        if (channel_type.lower() == "raw"):
           return RAW_CHANNELS[laser]
        elif (channel_type.lower() == "stitched"):
           return STITCHED_CHANNELS[laser]
        elif (channel_type.lower() == "corrected"):
           return CORRECTED_CHANNELS[laser]
        elif (channel_type.lower() == "downsized"):
           return DOWNSIZED_CHANNELS[laser]
        else:
           raise ValueError("Channel type must be raw, stitched, or corrected")

   def __check_sample_ind(self, sample_ind):
       if (sample_ind > len(self.get_samples())-1):
           raise ValueError("Only " + str(len(self.get_samples())-1) + " samples")  


