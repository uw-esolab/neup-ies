#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 09:23:58 2021

@author: gabrielsoto
"""

import pandas, os, json

class FileMethods(object):
    
    #defining directory where this particular file is located and other useful directories
    util_dir     = os.path.dirname(os.path.realpath(__file__))
    samsim_dir   = os.path.dirname(util_dir)
    neup_dir     = os.path.dirname(samsim_dir)
    parent_dir   = os.path.dirname(neup_dir)
    
    
    def read_csv_through_pandas(filepath):
        """ Method to read csv file and return data array
        
        Inputs:
            filepath (str) : full path to csv file
        Outputs:
            data_array (list) : data in either list (SSC_ARRAY) or nested list (SSC_MATRIX) form
        """
        dataframe = pandas.read_csv(filepath,header=None) #important: No Header assumed
        
        #recasting dataframe to a 1D numpy array if there is only 1 columnn in csv (SSC ARRAY)
        if dataframe.shape[1] == 1:
            data_array = dataframe.T.to_numpy()[0]
        #recasting dataframe to a ND numpy array, then to list of lists (equiv to SSC MATRIX)
        else:
            data_array = dataframe.to_numpy().tolist()
            
        return data_array


    def read_solar_resource_file(filepath, unit_registry):
        
        # setting a standard unit registry
        u = unit_registry 
        
        dataframe = pandas.read_csv(filepath,header=2)
        
        data = dataframe['Temperature']
        
        t_dry = data.to_numpy() * u.degC
        
        return t_dry.to('kelvin')


    def read_json(json_name):
        """ Method to read json file and return dictionaries
        
        Inputs:
            json_name (str) : name of json script found at 'neup-ies/simulations/json'
        Outputs:
            pysam_dict (dict) : dictionary of PySAM inputs + file names
            ssc_dict (dict) : dictionary of SSC inputs needed to run modules
            out_dict (dict) : dictionary of PySAM output keywords for extraction
        """
        #defining filepath for JSON script (FUTURE: add non-default functionality to input whole path)
        samsim_dir = FileMethods.samsim_dir
        json_filepath = os.path.join( samsim_dir, "json", json_name)
        json_filepath += ".json" 
        
        with open(json_filepath) as f:
            # loading JSON script to a dictionary
            D = json.load(f)
            
        #extracting specific dicts nested inside the JSON script
        pysam_dict = D['PySAM_inputs']
        ssc_dict   = D['SSC_inputs']
        out_dict   = D['PySAM_outputs']
        
        return pysam_dict, ssc_dict, out_dict

    def write_csv(data_list, columns, filename):
        """ Method to write csv file 
        
        Inputs:
            data_list (list of lists) : list of strings, each list having length n
            columns (list of str): list of column names with length n
            filename (str): name of file (without path)
            
        """

        params_dataframe = pandas.DataFrame(data_list,columns=columns)
        
        filepath = os.path.join(FileMethods.samsim_dir, "outputs", filename)
        
        params_dataframe.to_csv(filepath,index=False)

        