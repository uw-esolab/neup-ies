#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 11:37:50 2021

@author: gabrielsoto
"""

import PySAM.TcsmoltenSalt as TcsmoltenSalt
import modules.SolarTES as SolarTES
import PySAM.PySSC as pssc
import copy, time

# for some setup steps
soltes = SolarTES.SolarTES( is_dispatch=False )

# SSC dictionary
SSC_dict = soltes.SSC_dict

# SSC module name
plant_name = soltes.plant_name

# create plant data encoding for generic system
plant_datA = pssc.dict_to_ssc_table( SSC_dict, plant_name )

# =============================================================================
# Some helpful methods
# =============================================================================

def setUpPlant(Plant):
    # manually setting data arrays from csv files
    Plant.SolarResource.solar_resource_file         = soltes.solar_resource_file
    Plant.TimeOfDeliveryFactors.dispatch_factors_ts = soltes.df_array
    Plant.UserDefinedPowerCycle.ud_ind_od           = soltes.ud_array
    Plant.SystemControl.wlim_series                 = soltes.wl_array
    Plant.HeliostatField.helio_positions            = soltes.hp_array
    Plant.HeliostatField.eta_map                    = soltes.em_array
    Plant.HeliostatField.flux_maps                  = soltes.fm_array
    Plant.SystemControl.dispatch_series = [1.2]*8760

def safePlantExecute(Plant, name="Plant"):
    try: 
        Plant.execute()
        print("\n{0} successfully ran\n".format(name))
    except:
        print("\n{0} didn't run.\n".format(name))

def checkOutputDictionaryEmpty(dictionary, name="Plant"):
    print("{0} Output dictionary is {1}.".format( name, \
                            "NOT empty" if dictionary else "empty")    )
    
    return not bool(dictionary)

# =============================================================================
print("\n# ===========================================================")
print("# ====== Test #1: Copy Plant Data Encoding ======\n")
# =============================================================================

# copy of plant_dat
plant_datB = copy.deepcopy( plant_datA )

# create new Plant object
PlantA = TcsmoltenSalt.wrap(plant_datA)
outputA_pre = PlantA.Outputs.export()
print("--- Plant A created from PySSC plant data (originally from SSC dictionary).")

# new Plant object from copied plant_dat
PlantB = TcsmoltenSalt.wrap(plant_datB)
outputB_pre = PlantB.Outputs.export()
print("--- Plant B created from deep-copied PySSC Plant A data.")

# ================ execute Plant0 ================
setUpPlant(PlantA)
safePlantExecute(PlantA, "Plant A")
time.sleep(1)

# post-run Outputs of both plants
outputA_post = PlantA.Outputs.export()
outputB_post = PlantB.Outputs.export()

# ================ Checks ================
print("\n Conclusions: ")

# check if dictionaries are emtpy or not
ch = checkOutputDictionaryEmpty(outputA_pre,  "Plant A, pre-Plant A execution,")
ch = checkOutputDictionaryEmpty(outputA_post, "Plant A, post-Plant A execution,")
ch = checkOutputDictionaryEmpty(outputB_pre,  "Plant B, pre-Plant A execution,")
ch = checkOutputDictionaryEmpty(outputB_post, "Plant B, post-Plant A execution,")

# check that the execution of Plant worked
if outputA_pre != outputA_post:
    print("\nExecution of Plant A populated the Plant A Output dictionary.")
else:
    print("Execution of Plant A did NOT populate the Plant A Output dictionary.")

# check if PlantA and PlantB are still linked
if outputB_pre == outputB_post:
    print("Plant A and Plant B are NOT linked \n       (Plant B Output unaltered).")
else:
    print("Plant A and Plant B ARE linked \n       (Plant B Output was altered from Plant A execution).")


# =============================================================================
print("\n# ===========================================================")
print("# ====== Test #2: New Plant from Un-Executed Plant ======\n")
# =============================================================================

# create plant data encoding for generic system
plant_datC = pssc.dict_to_ssc_table( SSC_dict, plant_name )

# create new Plant object
PlantC = TcsmoltenSalt.wrap(plant_datC)
outputC_pre = PlantC.Outputs.export()
print("--- Plant C created from PySSC plant data (originally from SSC dictionary). Same as Plant A.")

# retrieve PlantC dictionary
plantC_dict_preExec = PlantC.export()

# create new Plant object from dictionary
PlantD = TcsmoltenSalt.new()
PlantD.assign(plantC_dict_preExec)
outputD_pre = PlantD.Outputs.export()
print("--- Plant D created from exported Plant C dictionary, before Plant C execution.")

# ================ execute Plant0 ================
setUpPlant(PlantC)
safePlantExecute(PlantC, "Plant C")
time.sleep(1)

# post-run Outputs of both plants
outputC_post = PlantC.Outputs.export()
outputD_post = PlantD.Outputs.export()

# ================ after execution ================

# retrieve PlantC dictionary
plantC_dict_postExec = PlantC.export()

# create new Plant object from dictionary
PlantE = TcsmoltenSalt.new()
PlantE.assign(plantC_dict_postExec)
outputE_post = PlantE.Outputs.export()
print("--- Plant E created from exported Plant C dictionary, after Plant C execution.")

# ================ Checks ================
print("\n Conclusions: ")

# check if dictionaries are emtpy or not
ch = checkOutputDictionaryEmpty(outputC_pre,  "Plant C, pre-Plant C execution,")
ch = checkOutputDictionaryEmpty(outputC_post, "Plant C, post-Plant C execution,")
ch = checkOutputDictionaryEmpty(outputD_pre,  "Plant D, pre-Plant C execution,")
ch = checkOutputDictionaryEmpty(outputD_post, "Plant D, post-Plant C execution,")
EisEmpty = checkOutputDictionaryEmpty(outputE_post, "Plant E, post-Plant C execution,")

# check that the execution of Plant worked
if outputC_pre != outputC_post:
    print("\nExecution of Plant C populated the Plant C Output dictionary.")
else:
    print("Execution of Plant C did NOT populate the Plant C Output dictionary.")

# check if PlantC and PlantD are still linked
if outputD_pre == outputD_post:
    print("Plant C and Plant D are NOT linked \n       (Plant D Output unaltered from Plant C execution).")
else:
    print("Plant C and Plant D ARE linked \n       (Plant D Output was altered from Plant C execution).")

if EisEmpty:
    print("Plant C and Plant E are NOT linked \n       (Plant E Output unaltered from Plant C execution).")
else:
    print("Plant C and Plant E ARE linked \n       (Plant E Output was altered from Plant C execution).")