#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:44:30 2021

@author: gabrielsoto
"""

import PySAM.NuclearTes as NuclearTes
from modules.GenericSSCModule import GenericSSCModule
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods

class NuclearTES(GenericSSCModule):
    
    def __init__(self, plant_name="nuclear_tes", json_name="model1"):
        
        # initialize Generic module, csv data arrays should be saved here
        GenericSSCModule.__init__( self, plant_name, json_name )


    def run_sim(self):
        """ Method to run single simulation for Generic System
        """
        
        # create Plant object and execute it
        plant = self.create_plant_object( )
        plant.execute( )
        
        # use executed Plant object to create Grid object and execute it
        grid  = self.create_grid_object( plant )
        grid.GridLimits.grid_curtailment = self.gc_array  #set curtailment to be really high
        grid.execute( )
        
        # use executed Plant object to create SingleOwner object and execute it
        so    = self.create_so_object( plant )
        so.execute( )
        
        return plant, grid, so
        
        
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
        
        
    def create_plant_object(self):
        """ Method to create Plant object
        """
        
        # create plant data encoding for generic system
        plant_dat = pssc.dict_to_ssc_table( self.SSC_dict, self.plant_name )
        
        # create new Plant object
        plant = NuclearTes.wrap(plant_dat)
        
        # manually setting data arrays from csv files
        plant.SolarResource.solar_resource_file         = self.solar_resource_file
        plant.TimeOfDeliveryFactors.dispatch_factors_ts = self.df_array
        plant.UserDefinedPowerCycle.ud_ind_od           = self.ud_array
        plant.SystemControl.wlim_series                 = self.wl_array
        plant.HeliostatField.helio_positions            = self.hp_array
        plant.HeliostatField.eta_map                    = self.em_array
        plant.HeliostatField.flux_maps                  = self.fm_array
        plant.SystemControl.dispatch_series = [1.2]*8760
        
        return plant
    
