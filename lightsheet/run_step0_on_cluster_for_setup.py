#DON'T ACTUALLY RUN THIS FILE, COPY AND PASTE LINES INTO TERMINAL
#lightsheet package
preprocessing.generateparamdict(os.getcwd(), **params) 
if not os.path.exists(os.path.join(params['outputdirectory'], 'lightsheet')): shutil.copytree(os.getcwd(), os.path.join(params['outputdirectory'], 'lightsheet'), ignore=shutil.ignore_patterns(*('.pyc','CVS','.git','tmp','.svn', 'TeraStitcher-Qt4-standalone-1.10.11-Linux'))) 


#clearmap package
updateparams(os.getcwd(), **params) # e.g. single job assuming directory_determiner function has been properly set
#copy folder into output for records
if not os.path.exists(os.path.join(params['outputdirectory'], 'clearmap_cluster')): shutil.copytree(os.getcwd(), os.path.join(params['outputdirectory'], 'clearmap_cluster'), ignore=shutil.ignore_patterns('^.git')) #copy run folder into output to save run info

