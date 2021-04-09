#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:40:22 2021

@author: gabrielsoto
"""

import PySAM.GenericSystem as GenericSystem
import PySAM.Grid as Grid
import PySAM.Singleowner as Singleowner
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods
import astropy.units as u
import numpy as np
import copy

class GenericSSCModule(object):
    
    def __init__(self, plant_name="generic_system", json_name="100mW_Generic"):
        
        # grab names, either default here or from child class
        self.json_name  = json_name
        self.plant_name = plant_name
        
        # read in dictionaries from json script
        PySAM_dict, SSC_dict = FileMethods.read_json( self.json_name )
        
        # storing SSC and Pyomo time horizons, inputs are in unit of hours
        self.ssc_horizon   = PySAM_dict['ssc_horizon'] * u.hr
        self.pyomo_horizon = PySAM_dict['pyomo_horizon'] * u.hr
        
        # save SSC_dict for usage later
        self.SSC_dict = SSC_dict
        
        # save csv arrays to class 
        self.store_csv_arrays( PySAM_dict )


    def run_sim(self):
        """ Method to run single simulation for Generic System
        """

        # create Plant object and execute it
        # TODO: here we can make a new method solely for plant execution
        #       - it could ask for a flag to run one time or in a loop?
        #       - OR we have a separate method that runs things in a loop
        #               - reason for this is that we need to log gen output
        #               - and pass it into grid and so
        self.create_Plant( )
        self.Plant.execute( )
        
        # use executed Plant object to create Grid object and execute it
        self.create_Grid( )
        self.Grid.execute( )
        
        # use executed Plant object to create SingleOwner object and execute it
        self.create_SO( )
        self.SO.execute( )

          
    def store_csv_arrays(self, input_dict):
        """ Method to get data from specified csv files and store in class
        
        Inputs:
            input_dict (dict) : dictionary with csv file names
        """
        
        # saving location of solar resource file for SSC input
        parent_dir = FileMethods.parent_dir
        self.solar_resource_file = parent_dir + input_dict['solar_resource_rel_parent']
        
        
    def create_Plant(self):
        """ Method to create Plant object for the first time
        """
        
        # create plant data encoding for generic system
        plant_dat = pssc.dict_to_ssc_table( self.SSC_dict, self.plant_name )
        
        # create new Plant object
        self.Plant = GenericSystem.wrap(plant_dat)

    
    def create_Grid(self):
        """ Method to create Grid object for the first time
        """
        
        # create grid data encoding for grid
        grid_dat = pssc.dict_to_ssc_table( self.SSC_dict, "grid" )
        
        # create new Grid object from existing Plant object
        self.Grid = Grid.from_existing( self.Plant )
        
        # import Grid-specific data to Grid object
        self.Grid.assign(Grid.wrap(grid_dat).export())


    def create_SO(self):
        """ Method to create SingleOwner object for the first time
        """
        
        # create singleowner data encoding for singleowner object
        so_dat   = pssc.dict_to_ssc_table( self.SSC_dict, "singleowner" )
        
        # create new Singleowner object from existing Plant object
        self.SO = Singleowner.from_existing( self.Plant )
        
        # import Singleowner-specific data to Singleowner object
        self.SO.assign(Singleowner.wrap(so_dat).export())

