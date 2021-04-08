#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:40:22 2021

@author: gabrielsoto
"""

from util.FileMethods import FileMethods

class GenericSSCModule(object):
    
    def __init__(self, plant_name="generic_system", json_name="100mW_Generic",
                       grid_name="grid", financial_name="singleownder" ):
        
        # grab json file name, either default here or from child class
        self.json_name  = json_name
        self.plant_name = plant_name
        self.grid_name  = grid_name
        self.financial_name = financial_name

        # read in dictionaries from json script
        PySAM_dict, SSC_dict = FileMethods.read_json( self.json_name )
        
        # save SSC_dict for usage later
        self.SSC_dict = SSC_dict
        
        # save csv arrays to class 
        self.get_csv_arrays( PySAM_dict )
        
        
    def get_csv_arrays(self, input_dict):
        """ Method to get data from specified csv files
        
        Inputs:
            input_dict (dict) : dictionary with csv file names
        """
        
        # saving location of solar resource file for SSC input
        parent_dir = FileMethods.parent_dir
        self.solar_resource_file = parent_dir + input_dict['solar_resource_rel_parent']
        
        
    def run_single_sim():
        
        
        return True