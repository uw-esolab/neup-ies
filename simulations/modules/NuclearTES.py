#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:44:30 2021

@author: gabrielsoto
"""

from modules.GenericSSCModule import GenericSSCModule
from util.FileMethods import FileMethods

class NuclearTES(GenericSSCModule):
    
    def __init__(self, plant_name="nuclear_tes", json_name="model1"):
        
        # initialize Generic module, csv data arrays should be saved here
        GenericSSCModule.__init__(self, plant_name, json_name)
        
        
    def get_csv_arrays(self, input_dict):
        """ Method to get data from specified csv files
        
        Inputs:
            input_dict (dict) : dictionary with csv file names
        """
        
        # saving location of solar resource file for SSC input
        parent_dir = FileMethods.parent_dir
        self.solar_resource_file = parent_dir + input_dict['solar_resource_rel_parent']
        
        # read csv and save data to arrays
        self.df_array = FileMethods.read_csv_through_pandas(input_dict['dispatch_factors_file'])
        self.ud_array = FileMethods.read_csv_through_pandas(input_dict['ud_file'])
        self.wl_array = FileMethods.read_csv_through_pandas(input_dict['wlim_file'])
        self.hp_array = FileMethods.read_csv_through_pandas(input_dict['helio_file'])
        self.gc_array = FileMethods.read_csv_through_pandas(input_dict['grid_file'])
        self.em_array = FileMethods.read_csv_through_pandas(input_dict['eta_file'])
        self.fm_array = FileMethods.read_csv_through_pandas(input_dict['flux_file'])
        
        
    def run_single_sim(self):
        return True
