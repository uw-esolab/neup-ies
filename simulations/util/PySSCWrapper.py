#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 16:00:31 2021

@author: gabrielsoto
"""

from util.FileMethods import FileMethods 
from scripts.PySSC_scripts.PySSC import PySSC
import os, json
import numpy as np
import pickle as pickle
import hashlib as hashlib

class PySSCWrapper(object):
    """
    The PySSCWrapper class is a util class used to run PySSC from a given JSON
    script (similar to what we're doing in the PySAM NE2 modules). The main 
    use-case for this method is mixed-mode debugging:
        (1) choose an existing script from the /scripts directory that calls this method
        (2) run the script
        (3) PySSCWrapper class is initialized
        (4) PySSCWrapper creates PySSC data structures and modules
        (5) PySSC calls SSC
        (6) We can debug SSC if we're running it through CodeLite or Eclipse
    
    For more information: https://github.com/uw-esolab/docs/blob/main/sam/debugSSCwithPySSC_Linux_CodeLiteIDE.md
    
    To ensure the debug version of SSC, look into the https://github.com/uw-esolab/neup-ies/build_debug_SAM bash script
    """
    
    def __init__(self, json_name='model1', is_debug=True):
        """ Initializes the PySSCWrapper module
        
        This method initialized the PySSCWrapper module. If we intend to run
        in debug mode, the flag used as input ensures we point to the correct
        version of the SSC library (libSSCd.so file which is the Debug version
        rather than the Release version). The method also finds and extracts
        information from the specified JSON script. 
        
        Inputs:
            json_name (str)  : name of JSON script with input data for module
            is_debug (bool)  : flag to determine if running in debug mode
        """
        
        self.is_debug = is_debug
        
        # defining useful directories
        self.samsim_dir = FileMethods.samsim_dir.encode("utf-8")
        self.NE2_dir    = FileMethods.parent_dir.encode("utf-8")
        
        # defining filepath for user selected JSON script
        self.json_filepath = os.path.join( FileMethods.samsim_dir, "json", json_name)
        self.json_filepath += ".json" 
        
        # loading JSON script to a dictionary
        with open(self.json_filepath) as f:
            D = json.load(f)
        
        #extracting specific dicts nested inside the JSON script
        self.pysam_dict = D['PySAM_inputs']
        self.sscdict    = D['SSC_inputs']


    def run_sim(self, run_dispatch_targets=False):
        """ Method to run full SSC simulation through PySSC
        
        This method creates an SSC API, saves input data from the specified
        JSON script to it and runs a full SSC simulation through PySSC. Much
        like PySAM, it creates a Plant, Grid, and Financial module and runs
        each in order. 
        """
        
        # creating SSC API and Data object(?) for PySSC
        self.create_API()
        
        # option to search for file with dispatch_targets and use arrays
        self.run_dispatch_targets = run_dispatch_targets
        
        #----Plant sims
        # initialize Plant-specific data
        plant_mod, plant_name = self.init_plant_module()
        # run Plant simulation
        plant_mod = self.run_module( plant_mod, plant_name )
        
        #----Grid sims
        # initialize Grid-specific data
        grid_mod, grid_name = self.init_grid_module()
        # run Grid simulation
        grid_mod = self.run_module( grid_mod, grid_name )
        
        #----SingleOwner sims
        # initialize SingleOwner-specific data
        fin_mod, fin_name = self.init_fin_module()
        # run SingleOwner simulation
        fin_mod = self.run_module( fin_mod, fin_name )

    
    def create_API(self):
        """ Method to create and save a unique SSC API to self
        
        This method also ensures we use the correct SSC library file, whether
        we need the Debug or Release version (both could be in play in different
        locations). 
        """
        
        self.sscapi = PySSC(self.is_debug)
        self.sscdata = self.sscapi.data_create()
    
    
    def init_plant_module(self):
        """ Method to initialize the Plant module within PySSC
        
        This method creates the Plant module using SSC data from csv files
        as well as from the SSC_dict within the specified JSON script. 

        Outputs:
            plant_module (int)       : some sort of data structure for the Plant (idk, PySSC stuff)
            plant_module_name (str)  : name of Plant module, or the cmod SSC file
        """
        
        slash = '/'.encode("utf-8")

        # file locations
        solar_resource_file   = self.NE2_dir    + slash + self.pysam_dict['solar_resource_rel_parent'].encode("utf-8")
        dispatch_factors_file = self.samsim_dir + slash + self.pysam_dict['dispatch_factors_file'].encode("utf-8")
        ud_file               = self.samsim_dir + slash + self.pysam_dict['ud_file'].encode("utf-8")
        wlim_file             = self.samsim_dir + slash + self.pysam_dict['wlim_file'].encode("utf-8")
        helio_file            = self.samsim_dir + slash + self.pysam_dict['helio_file'].encode("utf-8")
        eta_file              = self.samsim_dir + slash + self.pysam_dict['eta_file'].encode("utf-8")
        flux_file             = self.samsim_dir + slash + self.pysam_dict['flux_file'].encode("utf-8")
        
        # set data from files
        self.sscapi.data_set_string( self.sscdata, b'solar_resource_file', solar_resource_file )
        self.sscapi.data_set_matrix_from_csv( self.sscdata, b'eta_map',             eta_file)
        self.sscapi.data_set_matrix_from_csv( self.sscdata, b'flux_maps',           flux_file)
        self.sscapi.data_set_matrix_from_csv( self.sscdata, b'helio_positions',     helio_file)
        self.sscapi.data_set_matrix_from_csv( self.sscdata, b'ud_ind_od',           ud_file)
        self.sscapi.data_set_array_from_csv(  self.sscdata, b'wlim_series',         wlim_file)
        self.sscapi.data_set_array_from_csv(  self.sscdata, b'dispatch_factors_ts', dispatch_factors_file)
        
        # set data for dispatch targets if exists
        if self.run_dispatch_targets:
            hash_exists, hash_filepath = self.generate_hash()
            
            if hash_exists:
                print("\nDispatch Targets exists in {0}.\n".format(hash_filepath)) 
                
                # loading dispatch targets from file if it exists
                with open(hash_filepath, 'rb') as f:
                    self.DT = pickle.load(f)
                
                # setting dispatch targets to True
                self.sscapi.data_set_number(self.sscdata, 'is_dispatch_targets'.encode("utf-8"), 1)
                
                # setting dispatch target arrays in SSC inputs
                for key in self.DT.keys():
                    self.sscapi.data_set_array(self.sscdata, key.encode("utf-8"), self.DT[key])
                    
            else:
                dt = {'is_rec_su_allowed_in': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                      'is_rec_sb_allowed_in': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                      'is_pc_su_allowed_in': [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1], 
                      'is_pc_sb_allowed_in': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                      'q_pc_target_su_in': [0.0, 0.0, 357.5268818742051, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 357.5268818742051, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                      'q_pc_target_on_in': [0.0, 0.0, 874.4086, 874.4086, 874.4086, 874.4086, 1501.6129, 1501.6129, 1501.6129, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 592.47312, 1501.6129, 1501.6129, 1501.6129, 1501.6129, 1156.96818, 1156.69322, 1156.69322], 
                      'q_pc_max_in': [1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615, 1501.6129038716615], 
                      'f_nuc_to_tes_target': [1.0, 1.0, 0.07956989263157895, 0.07956989263157895, 0.07956989263157895, 0.07956989263157895, 0.0, 0.0, 0.0, 0.6905615321637427, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.37634408421052634, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}
                
                self.sscapi.data_set_number(self.sscdata, 'is_dispatch_targets'.encode("utf-8"), 1)
                
                for k in dt.keys():
                    tmp = np.hstack([dt[k], np.zeros(8760-24) ])
                    self.sscapi.data_set_array(self.sscdata, k.encode("utf-8"), tmp)
                
                print("\nDispatch Targets do not exist for this configuration. \n Will write to {0}.\n".format(hash_filepath))
            
        
        # setting dispatch series separately
        dispatch_series = [1.2]*8760
        self.sscapi.data_set_array( self.sscdata, b'dispatch_series',  dispatch_series)
        
        # cycle through JSON script, saving data as respective SSC Datatype, for Plant
        plant_module, plant_module_name = self.set_ssc_data_from_dict(0)
        
        return plant_module, plant_module_name


    def set_ssc_data_from_dict(self, start_ind):
        """ Method used to extract SSC data from a dictionary in a JSON script
        
        Method adapted from LORE team (particularly the loop). It checks the
        data type each entry of the JSON script SSC dictionary and uses the 
        appropriate SSC API method to store the data. 
        
        Inputs:
            start_ind (int)  : index of computer module (0: plant, 1: grid, 2: financial)
        Outputs:
            ssc_module (int)       : some sort of data structure for the Plant (idk, PySSC stuff)
            ssc_module_name (str)  : name of Plant module, or the cmod SSC file
        """
        
        # defining compute module name to extract data from
        module_name = 'compute_module_' + str(start_ind)
        
        # defining 'book-end' to end data extraction
        string_end  = 'compute_module_' + str(start_ind+1) if start_ind != 2 else 'number_compute_modules'
        
        # local instance of SSC objects
        sscdict = self.sscdict
        sscapi  = self.sscapi
        sscdata = self.sscdata
        
        # cycle through all SSC dicts from the JSON script
        for key in sscdict.keys():        
            
            # skip the first instance that holds the SSC module name
            if key == module_name:
                continue
            # end the loop if we hit the book-end string
            elif key == string_end:
                break
            # otherwise, try to get all values in between book-ends
            else:
                try: 
                    # set numbers if data is either int or float
                    if type(sscdict[key]) in [int, float]:
                        sscapi.data_set_number(sscdata, key.encode("utf-8"), sscdict[key])
                    # set boolean if data is boolean
                    elif type(sscdict[key]) == bool:
                        sscapi.data_set_number(sscdata, key.encode("utf-8"), 1 if sscdict[key] else 0)
                    # set string if data is a str
                    elif type(sscdict[key]) == str:
                        sscapi.data_set_string(sscdata, key.encode("utf-8"), sscdict[key].encode("utf-8"))
                    # set matrix or array if data is a list
                    elif type(sscdict[key]) == list:
                        if len(sscdict[key]) > 0:
                            if type(sscdict[key][0]) == list:
                                sscapi.data_set_matrix(sscdata, key.encode("utf-8"), sscdict[key])
                            else:
                                sscapi.data_set_array(sscdata, key.encode("utf-8"), sscdict[key])
                        else:
                            #print ("Did not assign empty array " + key)
                            pass           
                    # re-run loop for this key if data entry is a dict
                    elif type(sscdict[key]) == dict:
                        # TODO: fix this, so far there are no dicts defined within the SSC dict
                        table = sscapi.data_create()
                        self.set_sscdata_from_dict(sscapi, table, sscdict[key])
                        sscapi.data_set_table(sscdata, key.encode("utf-8"), table)
                        sscapi.data_free(table)
                    # raise error if key not found
                    else:
                       print ("Could not assign variable " + key )
                       raise KeyError
                # raise error for unknown data type     
                except:
                    print ("Error assigning variable " + key + ": bad data type")
        
        # save module name and module object to return
        ssc_module_name = sscdict[module_name]
        ssc_module = sscapi.module_create(ssc_module_name.encode("utf-8"))
        
        # update SSC objects saved to self
        self.sscdict = sscdict
        self.sscapi  = sscapi
        self.sscdata = sscdata
        
        return ssc_module, ssc_module_name
        
        
    def run_module(self, module, name):
        """ Method used to run a given SSC module through PySSC
        
        Method adapted from LORE team (particularly the loop). It runs the given
        module (Plant, Grid, Financial) by calling SSC through the SSC API. 
        
        Inputs:
            module (int)       : some sort of data structure for the Plant (idk, PySSC stuff)
            module_name (str)  : name of Plant module, or the cmod SSC file
        Outputs:
            module (int)       : updated module with SSC results
        """
        
        sscdata = self.sscdata
        sscapi  = self.sscapi
        
        sscapi.module_exec_set_print( 0 )
        if sscapi.module_exec(module, sscdata) == 0:
            print ('    ' + name + ' simulation error')
            idx = 1
            msg = sscapi.module_log(module, 0)
            while (msg != None):
                print ('    : ' + msg.decode("utf - 8"))
                msg = sscapi.module_log(module, idx)
                idx = idx + 1
            SystemExit( "Simulation Error" )
        sscapi.module_free(module)
        
        self.sscapi  = sscapi
        self.sscdata = sscdata
        
        return module
        
    
    def init_grid_module(self):
        """ Method to initialize the Grid module within PySSC
        
        This method creates the Grid module using SSC data from csv files
        as well as from the SSC_dict within the specified JSON script. 

        Outputs:
            grid_module (int)       : some sort of data structure for the Grid mod (idk, PySSC stuff)
            grid_module_name (str)  : name of Grid module, or the cmod SSC file
        """
        
        slash = '/'.encode("utf-8")
        
        # file location
        grid_file  = self.samsim_dir + slash + self.pysam_dict['grid_file'].encode("utf-8")
        
        # get array from csv file for grid curtailment
        self.sscapi.data_set_array_from_csv(  self.sscdata, b'grid_curtailment', grid_file)
        
        # cycle through JSON script, saving data as respective SSC Datatype, for Grid
        grid_mod, grid_name = self.set_ssc_data_from_dict(1)
        
        return grid_mod, grid_name


    def init_fin_module(self):
        """ Method to initialize the Financial module within PySSC
        
        This method creates the Financial module using SSC data from the 
        SSC_dict within the specified JSON script. 

        Outputs:
            fin_module (int)       : some sort of data structure for the Financial mod (idk, PySSC stuff)
            fin_module_name (str)  : name of Financial module, or the cmod SSC file
        """
        
        # cycle through JSON script, saving data as respective SSC Datatype, for Grid
        fin_mod, fin_name = self.set_ssc_data_from_dict(2)
        
        return fin_mod, fin_name 


    def generate_hash(self):
        """ Method to create unique hash for given JSON inputs
        
        This method creates a unique, permanent hash for a given JSON script.
        That is, it gathers all of the JSON inputs (including SSC and PySAM inputs)
        from the designated script and converts both their keynames and values
        to strings. It collects all of these into a single string variable and
        then creates a new hexadecimal string or "hash" for that giant string. 
        This serves as a unique identifier or "fingerprint" for all the values 
        in the JSON script. This is then used later on as the file name containing
        outputs from this particular run. Any small changes to the JSON script 
        will result in a drastically different hash, and therefore a new output file.
        If a simulation has already been run with the given JSON script, it can
        just pull results from the already created hash file instead of needlessly
        repeating the simulation. 

        Returns:
            hash_exists (bool): 
                if True, a hash file currently exists with all given JSON inputs
            filepath (str): 
                absolute filepath to the hash file in outputs directory
        """
        
        # static start of the filename
        filename = "NuclearTES__"   if self.sscdict['compute_module_0'] == 'nuclear_tes' else \
                   "SolarTES__"     if self.sscdict['compute_module_0'] == 'tcsmolten_salt' else \
                   "DualPlantTES__" if self.sscdict['compute_module_0'] == 'nuclear_mspt_tes' else \
                   "__"
        
        if self.sscdict['compute_module_0'] == 'nuclear_mspt_indirect_tes':
            filename = 'IndirectNuclearTES__' if self.sscdict['q_dot_rec_des'] == 0 else 'DualIndirectTES__'
        
        # initializing empty string
        extstr = ''
        
        # adding SSC dictionary names and values to the existing string
        sscdict = self.sscdict
        for s in sscdict.keys():
            extstr += "{0}: {1} ".format( s, str(sscdict[s]) )
        
        # adding PySAM dictionary names and values to the existing string
        pysamdict = self.pysam_dict
        for p in pysamdict.keys():
            extstr += "{0}: {1} ".format( p, str(pysamdict[p]) )
        
        # creating a unique hash from giant string of values using md5 algorithm
        json_hash = hashlib.md5( extstr.encode('utf-8') ).hexdigest()
        
        # adding the unique hexadecimal hash to the starting filename string
        filename += json_hash
        
        # creating full path to output hash file
        filepath = os.path.join( FileMethods.output_dir , filename + '.dispatchTargets')
        
        # checking if this current hash exists already
        hash_exists = os.path.exists(filepath)
        
        return hash_exists, filepath
    
    
    def get_array(self, out_str):
        """ Method to extract data from SSC results
        
        This method gets SSC outputs and converts it to arrays.

        Inputs:
            out_str (str)       : string name of the desired output
        Outputs:
            out_array (ndarray) : SSC result in array form
        """
        ssc = self.sscapi
        data = self.sscdata
        
        out_array = ssc.data_get_array(data,out_str.encode('utf-8'))
        out_array = np.asarray(out_array)
        return out_array
        

    
if __name__ == "__main__": 
    pw = PySSCWrapper()
    pw.run_sim()