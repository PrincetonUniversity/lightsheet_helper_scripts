#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 10:36:05 2019

@author: wanglab
"""

import pickle

with open("/home/wanglab/Documents/cfos_inputs/cfos_points_dictionary.p", "rb") as f:
    w = pickle.load(f, encoding = "latin1")

pickle.dump(w, open("/home/wanglab/Documents/cfos_inputs/cfos_points_py2_dictionary.p","wb"), protocol=2)