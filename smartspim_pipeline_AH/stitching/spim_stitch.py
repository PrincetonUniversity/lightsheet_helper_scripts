import os,sys
from subprocess import run,PIPE


def stitch_step0(volin):
    """ call function for terastitcher """
    cmd = "terastitcher -1 --volin={} --ref1=x --ref2=y --ref3=z --vxl1=1.81 --vxl2=1.81 --vxl3=2 --projout=xml_import".format(
        volin)
    print("running command: ",cmd)
    result = run(cmd,
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
    """ terastitcher uses stdout for stderr and doesn't return a proper exit code
    so to catch errors we just need to check if stdout has anything in it """
    print(result)
    if b'ERROR' in result.stdout:
        raise Exception(result.stdout)
   
    return 

def stitch_step1(volin):
    """ call function for terastitcher """
    projin = os.path.join(volin, "xml_import.xml")
    projout = os.path.join(volin,"xml_displcomp.xml")
    cmd = "terastitcher --displcompute --projin={} --subvoldim=100 --projout={}".format(
        projin,projout) 
    print("running command: ",cmd)
    
    result = run(cmd,
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
    """ terastitcher uses stdout for stderr and doesn't return a proper exit code
    so to catch errors we just need to check if stdout has anything in it """
    if b'ERROR' in result.stdout:
        raise Exception(result.stdout)
    
    return 

def stitch_step2(volin):
    """ Runs the projection steps """
    
    """ Part (1/3) """
    projin1 = os.path.join(volin, "xml_displcomp.xml")
    cmd1 = "terastitcher --displproj --projin={}".format(
        projin1) 
    print("running command: ",cmd1)

    result1 = run(cmd1,
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
   
    if b'ERROR' in result1.stdout:
        raise Exception("error found in terastitcher output: {}".format(result1.stdout))
    
    """ Part (2/3) """
    projin2 = os.path.join(volin, "xml_displproj.xml")
    projout2 = os.path.join(volin,"xml_displthres.xml")
    cmd2 = "terastitcher --displthres --projin={} --projout={} --threshold=0.7".format(
        projin2,projout2) 
    print("running command: ",cmd2)

    result2 = run(cmd2,
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
   
    if b'ERROR' in result2.stdout:
        raise Exception("error found in terastitcher output: {}".format(result2.stdout))   

    """ Part (3/3) """
    projin3 = os.path.join(volin,"xml_displthres.xml")
    projout3 = os.path.join(volin,"xml_placetiles.xml")
    cmd3 = "terastitcher --placetiles --projin={} --projout={} --algorithm=MST".format(
        projin3,projout3) 
    print("running command: ",cmd3)

    result3 = run(cmd3,
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
   
    if b'ERROR' in result3.stdout:
        raise Exception("error found in terastitcher output: {}".format(result3.stdout))   

    # sp_call("terastitcher --merge --projin=%s --volout=%s --imout_depth=16 --resolutions=0" % (os.path.join(volin, "xml_placetiles.xml"), dest))
    
    return 

def stitch_step3(volin,volout):
    """ Runs the projection steps """
    
    """ Part (1/3) """
    projin = os.path.join(volin, "xml_placetiles.xml")
    cmd = "terastitcher --merge --projin={} --volout={} --imout_depth=16 --resolutions=0".format(
        projin,volout) 
    print("running command: ",cmd)

    result = run(cmd,
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
   
    if b'ERROR' in result.stdout:
        raise Exception("error found in terastitcher output: {}".format(result.stdout))
    
    return 

if __name__ == "__main__":

    #get jobids from SLURM or argv
    print(sys.argv)
    stepid = sys.argv[1]
    input_dir = sys.argv[2]
    output_dir = sys.argv[3]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    #######################STEP 0 #######################
    #Make parameter dictionary and setup destination
    #####################################################
    if stepid == "step0":
        ###make parameter dictionary and pickle file:
        print("running step 0!")
        stitch_step0(volin=input_dir)
    elif stepid == "step1":
        ###make parameter dictionary and pickle file:
        print("running step 1!")
        stitch_step1(volin=input_dir)
    elif stepid == "step2":
        ###make parameter dictionary and pickle file:
        print("running step 2!")
        stitch_step2(volin=input_dir)
    elif stepid == "step3":
        print("running step 3!")
        stitch_step3(volin=input_dir,volout=output_dir)
