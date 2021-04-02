#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 09:23:58 2021

@author: gabrielsoto
"""

import pandas, os

class FileMethods(object):
    
    #defining directory where this particular file is located and other useful directories
    util_dir      = os.path.dirname(os.path.realpath(__file__))
    neup_dir      = os.path.dirname(util_dir)
    scripts_dir   = os.path.join(neup_dir,"python-scripts")
    test_dir      = os.path.join(neup_dir,"tests")
    parent_dir    = os.path.dirname(neup_dir)
    json_test_dir = os.path.join(scripts_dir,"json-scripts/tests")
    
    def read_csv_through_pandas(filepath):
        """ Method to read csv file and return data array
        
        Inputs:
            filepath (str) - full path to csv file
        Outputs:
            data_array (list) - data in either list (SSC_ARRAY) or nested list (SSC_MATRIX)
        
        """
        dataframe = pandas.read_csv(filepath,header=None) #important: No Header assumed
        #recasting dataframe to a 1D numpy array if there is only 1 columnn in csv (SSC ARRAY)
        if dataframe.shape[1] == 1:
            data_array = dataframe.T.to_numpy()[0]
        #recasting dataframe to a ND numpy array, then to list of lists (equiv to SSC MATRIX)
        else:
            data_array = dataframe.to_numpy().tolist()
        return data_array