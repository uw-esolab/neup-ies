#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:47:58 2020

Most recently tested against PySAM 2.1.4

@author: frohro
"""
import json, os
import pandas
import PySAM.NuclearTes as NuclearTes
import PySAM.Grid as Grid
import PySAM.Singleowner as Singleowner
import PySAM.PySSC as pssc
ssc = pssc.PySSC()

# defining directories
cwd        = os.getcwd() #should we make this static? cwd can change based on IDE preferences
neup_dir   = os.path.dirname(cwd)
parent_dir = os.path.dirname(neup_dir)
ssc_dir    = parent_dir + '/build_ssc/ssc/libssc.so'

# method to read in csv file and output data array in correct shape
def read_pandas(filepath):
    dataframe = pandas.read_csv(filepath,header=None)
    if dataframe.shape[1] == 1:
        data_array = dataframe.T.to_numpy()[0]
    else:
        data_array = dataframe.to_numpy().tolist()
    return data_array

# solar file that comes with SAM repository
solar_resource_file = parent_dir + '/sam/deploy/solar_resource/tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv'

# creating data arrays from csv files
df_array = read_pandas(cwd + '/data-files/dispatch_factors_ts.csv')
ud_array = read_pandas(cwd + '/data-files/ud_ind_od.csv')
wl_array = read_pandas(cwd + '/data-files/wlim_series.csv')
hp_array = read_pandas(cwd + '/data-files/helio_positions.csv')
gc_array = read_pandas(cwd + '/data-files/grid_curtailment.csv')
em_array = read_pandas(cwd + '/data-files/eta_map.csv')
fm_array = read_pandas(cwd + '/data-files/flux_maps.csv')

# defining modules to run
with open("json-scripts/model1.json") as f:
    # loading json script to a dictionary
    dic = json.load(f)
    nt_dat = pssc.dict_to_ssc_table(dic, "nuclear_tes")
    grid_dat = pssc.dict_to_ssc_table(dic, "grid")
    so_dat = pssc.dict_to_ssc_table(dic, "singleowner")
    
    # creating Nuclear module from data
    nt = NuclearTes.wrap(nt_dat)

    # manually setting data arrays from csv files
    nt.SolarResource.solar_resource_file         = solar_resource_file
    nt.TimeOfDeliveryFactors.dispatch_factors_ts = df_array
    nt.UserDefinedPowerCycle.ud_ind_od           = ud_array
    nt.SystemControl.wlim_series                 = wl_array
    nt.HeliostatField.helio_positions            = hp_array
    nt.HeliostatField.eta_map                    = em_array
    nt.HeliostatField.flux_maps                  = fm_array
    nt.SystemControl.dispatch_series = [1.2]*8760
    
    # creating Grid module from existing Nuclear module
    grid = Grid.from_existing(nt)
    grid.assign(Grid.wrap(grid_dat).export())
    grid.GridLimits.grid_curtailment = gc_array #setting grid curtailment data taken from csv file

    # to create GenericSystem and Singleowner combined simulation, sharing the same data
    so = Singleowner.from_existing(nt)
    so.assign(Singleowner.wrap(so_dat).export())


nt.execute()
grid.execute()
so.execute()
print('Made it past execute.')
#print(nt.Outputs.export())  # as dictionary

annual_energy          = nt.Outputs.annual_energy/1e9
capacity_factor        = nt.Outputs.capacity_factor
annual_total_water_use = nt.Outputs.annual_total_water_use
ppa                    = so.Outputs.ppa
lppa_nom               = so.Outputs.lppa_nom
lppa_real              = so.Outputs.lppa_real
lcoe_nom               = so.Outputs.lcoe_nom
lcoe_real              = so.Outputs.lcoe_real
project_return_aftertax_npv = so.Outputs.project_return_aftertax_npv/1e6
flip_actual_irr        = so.Outputs.flip_actual_irr
flip_actual_year       = so.Outputs.flip_actual_year
project_return_aftertax_irr = so.Outputs.project_return_aftertax_irr
cost_installed         = so.Outputs.cost_installed/1e6
size_of_equity         = so.Outputs.size_of_equity/1e6
size_of_debt           = so.Outputs.size_of_debt/1e6

print('')
print('                        Nuclear')
print ('Annual energy (year 1)        =   ', annual_energy, ' TWh')
print ('Capacity factor (year 1)      =   ', capacity_factor, ' %')
print ('Annual Water Usage            =   ', annual_total_water_use, ' m3')
print ('PPA price (year 1)            =   ', ppa, ' ¢/kWh')
print ('Levelized PPA price (nominal) =   ', lppa_nom, ' ¢/kWh')
print ('Levelized PPA price (real)    =   ', lppa_real, ' ¢/kWh')
print ('Levelized COE (nominal)       =   ', lcoe_nom, ' ¢/kWh')
print ('Levelized COE (real)          =   ', lcoe_real, ' ¢/kWh')
print ('Net present value             =  $', project_return_aftertax_npv, ' M')
print ('Internal rate of return (IRR) =   ', flip_actual_irr, ' %')
print ('Year IRR is achieved          =   ', flip_actual_year)
print ('IRR at end of project         =   ', project_return_aftertax_irr, ' %')
print ('Net capital cost              =  $', cost_installed, ' M')
print ('Equity                        =  $', size_of_equity, ' M')
print ('Size of debt                  =  $', size_of_debt, ' M')
