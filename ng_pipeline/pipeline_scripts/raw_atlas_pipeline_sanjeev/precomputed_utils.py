from collections import defaultdict


def calculate_chunks(downsample, mip):
    """
    Chunks default to 64,64,64 so we want different chunks at 
    different resolutions
    """
    d = defaultdict(dict)
    d['full'][-1] = [1024,1024,1]
    d['full'][0] = [128,128,64]
    d['full'][1] = [128,128,64]
    d['full'][2] = [128,128,64]
    d['full'][3] = [128,128,64]
    d['full'][4] = [128,128,64]
    d['full'][5] = [64,64,64]
    d['full'][6] = [64,64,64]
    d['full'][7] = [64,64,64]
    d['full'][8] = [64,64,64]
    d['full'][9] = [64,64,64]

    try:
        result = d[downsample][mip]
    except:
        result = [64,64,64]
    return result

def calculate_factors(downsample, mip):
    """
    Scales get calculated by default by 2x2x1 downsampling
    """
    d = defaultdict(dict)
    d['full'][0] = [2,2,1]
    d['full'][1] = [2,2,2]
    d['full'][2] = [2,2,2]
    d['full'][3] = [2,2,2]
    d['full'][4] = [2,2,2]
    d['full'][5] = [2,2,2]
    d['full'][6] = [2,2,2]
    d['full'][7] = [2,2,2]
    d['full'][8] = [2,2,2]
    d['full'][9] = [2,2,2]

    try:
        result = d[downsample][mip]
    except:
        result = [2,2,1]
    return result
