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

cwd        = os.getcwd()
neup_dir   = os.path.dirname(cwd)
parent_dir = os.path.dirname(neup_dir)
ssc_dir    = parent_dir + '/build_ssc/ssc/libssc.so'

def read_pandas(filepath):
    dataframe = pandas.read_csv(filepath,header=None)
    if dataframe.shape[1] == 1:
        data_array = dataframe.T.to_numpy()[0]
    else:
        data_array = dataframe.to_numpy().tolist()
    return data_array

solar_resource_file = parent_dir + '/sam/deploy/solar_resource/tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv'

df_array = read_pandas(cwd + '/data-files/dispatch_factors_ts.csv')
ud_array = read_pandas(cwd + '/data-files/ud_ind_od.csv')
wl_array = read_pandas(cwd + '/data-files/wlim_series.csv')
hp_array = read_pandas(cwd + '/data-files/helio_positions.csv')
gc_array = read_pandas(cwd + '/data-files/grid_curtailment.csv')
em_array = read_pandas(cwd + '/data-files/eta_map.csv')
fm_array = read_pandas(cwd + '/data-files/flux_maps.csv')

with open("json-scripts/model1.json") as f:
    dic = json.load(f)
    gs_dat = pssc.dict_to_ssc_table(dic, "nuclear_tes")
    grid_dat = pssc.dict_to_ssc_table(dic, "grid")
    so_dat = pssc.dict_to_ssc_table(dic, "singleowner")

    gs = NuclearTes.wrap(gs_dat)

    # manually setting file locations
    gs.SolarResource.solar_resource_file         = solar_resource_file
    gs.TimeOfDeliveryFactors.dispatch_factors_ts = df_array
    gs.UserDefinedPowerCycle.ud_ind_od           = ud_array
    gs.SystemControl.wlim_series                 = wl_array
    gs.HeliostatField.helio_positions            = hp_array
    gs.HeliostatField.eta_map                    = em_array
    gs.HeliostatField.flux_maps                  = fm_array
    gs.SystemControl.dispatch_series = [1.2]*8760

    grid = Grid.from_existing(gs)
    grid.assign(Grid.wrap(grid_dat).export())

    # to create GenericSystem and Singleowner combined simulation, sharing the same data
    so = Singleowner.from_existing(gs)
    so.assign(Singleowner.wrap(so_dat).export())


gs.execute()
grid.execute()
so.execute()
print('Made it past execute.')
#print(gs.Outputs.export())  # as dictionary
