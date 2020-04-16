#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 18:43:25 2019

@author: wanglab
"""

from subprocess import check_output

#function to run
def sp_call(call):
    """ command line call function """ 
    print(check_output(call, shell=True)) 
    return


fls = ['/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk06/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tp02/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tp01/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tp07/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tp06/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk01/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tp08/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk05/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk07/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tpalpha/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk03/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tp09/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk11/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk02/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk04/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk08/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk10/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*',
 '/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_tpbeta/clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters*']

for fl in fls:
    call = "chmod 777 {}".format(fl)
    print(fl)
    try:
        sp_call(call)
        print(call)
    except:
        print("run registration")