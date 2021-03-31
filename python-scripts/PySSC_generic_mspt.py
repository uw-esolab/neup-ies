"""
Python script written for the NE-2 project.
Written to run simulations of SSC code through PySSC. This script is based on
    the auto-generated SAM script for Python3. 
Contains plotting utilities and is capable of running new models such as
    the Nuclear + TES model for the NE-2 project.

Authors: Gabriel Soto
"""
from PySSC import PySSC
import matplotlib.pyplot as plt
import numpy as np
import os

from pylab import *
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)

if __name__ == "__main__":
    ssc = PySSC()
    
    # Assumes that the SAM build folder is in the parent directory that houses neup-ies repo
    #    /parent_dir/build/..
    #    /parent_dir/neup-ies/python-scripts
    pdir = ssc.parent_dir 
    cwd  = pdir + '/neup-ies/python-scripts'
    print ('Current folder = ', cwd)
    print ('SSC Version = ', ssc.version())
    print ('SSC Build Information = ', ssc.build_info().decode("utf - 8"))
    ssc.module_exec_set_print(0)
    # getting process ID to attach in a C++ debugger
    pid = os.getpid()
    print("PID = ", pid)
    
    # File locations
    bprnt_dir = pdir.encode("utf-8") #parent directory in unicode coding
    bdata_dir = cwd.encode("utf-8")  #data repo directory in unicode conding
    solar_resource_file = bprnt_dir + b'/sam/deploy/solar_resource/tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv'
    dispatch_factors_ts = bdata_dir + b'/data-files/dispatch_factors_ts.csv'
    ud_ind_od           = bdata_dir + b'/data-files/ud_ind_od.csv'
    wlim_series         = bdata_dir + b'/data-files/wlim_series.csv'
    helio_positions     = bdata_dir + b'/data-files/helio_positions.csv'
    grid_curtailment    = bdata_dir + b'/data-files/grid_curtailment.csv'
    eta_map             = bdata_dir + b'/data-files/eta_map.csv'
    flux_maps           = bdata_dir + b'/data-files/flux_maps.csv'
    
    # =============================================================================
    # Starting data collection for run
    # =============================================================================
    data = ssc.data_create()

    # solarpilot ----------------------------------------------------------
    ssc.data_set_number( data, b'field_model_type', 2 )
    #-------------------------------------------------------------------------------
    
    # choosing nuclear vs csp ------------------------------------------------------
    module_name = b'tcsmolten_salt'
    #--------------------------------------------------------------------------------
    
    # changing end time of simulation -----------------------------------------------
    hours        = 24*7 # hours of simulation desired
    is_full_year = 1  # set to 1 for true if desired simulation time is 1 year
    if is_full_year:
        hours = int(31536000 / 3600) # 1 full year in hours
    ssc.data_set_number( data, b'time_start', 0 )
    ssc.data_set_number( data, b'time_stop', 3600*hours )
    #--------------------------------------------------------------------------------
    
    # setting dispatch parameters ---------------------------------------------------
    ssc.data_set_number( data, b'is_dispatch', 0) #set to 1 == True
    #--------------------------------------------------------------------------------
    
    ssc.data_set_string( data, b'solar_resource_file', solar_resource_file )
    ssc.data_set_number( data, b'ppa_multiplier_model', 0 )
    ssc.data_set_array_from_csv( data, b'dispatch_factors_ts', dispatch_factors_ts)
    ssc.data_set_number( data, b'gross_net_conversion_factor', 0.90000000000000002 )
    ssc.data_set_number( data, b'helio_width',              12.199999999999999 )
    ssc.data_set_number( data, b'helio_height',             12.199999999999999 )
    ssc.data_set_number( data, b'helio_optical_error_mrad', 1.53 )
    ssc.data_set_number( data, b'helio_active_fraction',    0.98999999999999999 )
    ssc.data_set_number( data, b'dens_mirror',              0.96999999999999997 )
    ssc.data_set_number( data, b'helio_reflectance',        0.90000000000000002 )
    ssc.data_set_number( data, b'rec_absorptance',          0.93999999999999995 )
    ssc.data_set_number( data, b'rec_hl_perm2',             30 )
    ssc.data_set_number( data, b'land_max', 9.5 )
    ssc.data_set_number( data, b'land_min', 0.75 )
    ssc.data_set_number( data, b'dni_des', 950 )
    ssc.data_set_number( data, b'p_start',                  0.025000000000000001 )
    ssc.data_set_number( data, b'p_track',                  0.055 )
    ssc.data_set_number( data, b'hel_stow_deploy',          8 )
    ssc.data_set_number( data, b'v_wind_max', 15 )
    ssc.data_set_number( data, b'c_atm_0', 0.0067889999999999999 )
    ssc.data_set_number( data, b'c_atm_1', 0.1046 )
    ssc.data_set_number( data, b'c_atm_2', -0.017000000000000001 )
    ssc.data_set_number( data, b'c_atm_3', 0.0028449999999999999 )
    ssc.data_set_number( data, b'n_facet_x', 2 )
    ssc.data_set_number( data, b'n_facet_y', 8 )
    ssc.data_set_number( data, b'focus_type', 1 )
    ssc.data_set_number( data, b'cant_type', 1 )
    ssc.data_set_number( data, b'n_flux_days', 8 )
    ssc.data_set_number( data, b'delta_flux_hrs', 2 )
    ssc.data_set_number( data, b'water_usage_per_wash',     0.69999999999999996 )
    ssc.data_set_number( data, b'washing_frequency', 63 )
    ssc.data_set_number( data, b'check_max_flux', 0 )
    ssc.data_set_number( data, b'sf_excess', 1 )
    ssc.data_set_number( data, b'tower_fixed_cost', 3000000 )
    ssc.data_set_number( data, b'tower_exp', 0.011299999999999999 )
    ssc.data_set_number( data, b'rec_ref_cost', 103000000 )
    ssc.data_set_number( data, b'rec_ref_area', 1571 )
    ssc.data_set_number( data, b'rec_cost_exp', 0.69999999999999996 )
    ssc.data_set_number( data, b'site_spec_cost', 16 )
    ssc.data_set_number( data, b'heliostat_spec_cost', 140 )
    ssc.data_set_number( data, b'plant_spec_cost', 1040 )
    ssc.data_set_number( data, b'bop_spec_cost', 290 )
    ssc.data_set_number( data, b'tes_spec_cost', 22 )
    ssc.data_set_number( data, b'land_spec_cost', 10000 )
    ssc.data_set_number( data, b'contingency_rate', 7 )
    ssc.data_set_number( data, b'sales_tax_rate', 5 )
    ssc.data_set_number( data, b'sales_tax_frac', 80 )
    ssc.data_set_number( data, b'cost_sf_fixed', 0 )
    ssc.data_set_number( data, b'fossil_spec_cost', 0 )
    ssc.data_set_number( data, b'flux_max', 1000 )
    ssc.data_set_number( data, b'opt_init_step', 0.059999999999999998 )
    ssc.data_set_number( data, b'opt_max_iter', 200 )
    ssc.data_set_number( data, b'opt_conv_tol', 0.001 )
    ssc.data_set_number( data, b'opt_flux_penalty', 0.25 )
    ssc.data_set_number( data, b'opt_algorithm', 1 )
    ssc.data_set_number( data, b'csp.pt.cost.epc.per_acre', 0 )
    ssc.data_set_number( data, b'csp.pt.cost.epc.percent', 13 )
    ssc.data_set_number( data, b'csp.pt.cost.epc.per_watt', 0 )
    ssc.data_set_number( data, b'csp.pt.cost.epc.fixed', 0 )
    ssc.data_set_number( data, b'csp.pt.cost.plm.percent', 0 )
    ssc.data_set_number( data, b'csp.pt.cost.plm.per_watt', 0 )
    ssc.data_set_number( data, b'csp.pt.cost.plm.fixed', 0 )
    ssc.data_set_number( data, b'csp.pt.sf.fixed_land_area', 45 )
    ssc.data_set_number( data, b'csp.pt.sf.land_overhead_factor', 1 )
    ssc.data_set_number( data, b'T_htf_cold_des', 290 )
    ssc.data_set_number( data, b'T_htf_hot_des', 574 )
    ssc.data_set_number( data, b'P_ref',                     115 )
    ssc.data_set_number( data, b'design_eff',                0.41199999999999998 )
    ssc.data_set_number( data, b'tshours', 10 )
    ssc.data_set_number( data, b'solarm', 2.3999999999999999 )
    ssc.data_set_number( data, b'N_panels', 20 )
    ssc.data_set_number( data, b'd_tube_out', 40 )
    ssc.data_set_number( data, b'th_tube', 1.25 )
    ssc.data_set_number( data, b'mat_tube', 2 )
    ssc.data_set_number( data, b'rec_htf', 17 )
    field_fl_props = [[ 0,   0,   0,   0,   0,   0,   0 ]]
    ssc.data_set_matrix( data, b'field_fl_props', field_fl_props )
    ssc.data_set_number( data, b'Flow_type', 1 )
    ssc.data_set_number( data, b'epsilon', 0.88 )
    ssc.data_set_number( data, b'hl_ffact', 1 )
    ssc.data_set_number( data, b'f_rec_min', 0.25 )
    ssc.data_set_number( data, b'rec_su_delay', 0.20000000000000001 )
    ssc.data_set_number( data, b'rec_qf_delay', 0.25 )
    ssc.data_set_number( data, b'csp.pt.rec.max_oper_frac', 1.2 )
    ssc.data_set_number( data, b'eta_pump', 0.84999999999999998 )
    ssc.data_set_number( data, b'piping_loss', 10200 )
    ssc.data_set_number( data, b'piping_length_mult', 2.6000000000000001 )
    ssc.data_set_number( data, b'piping_length_const', 0 )
    ssc.data_set_number( data, b'is_rec_model_trans', 0 )
    ssc.data_set_number( data, b'is_rec_startup_trans', 0 )
    ssc.data_set_number( data, b'rec_tm_mult', 1 )
    ssc.data_set_number( data, b'riser_tm_mult', 1 )
    ssc.data_set_number( data, b'downc_tm_mult', 1 )
    ssc.data_set_number( data, b'u_riser', 4 )
    ssc.data_set_number( data, b'th_riser', 15 )
    ssc.data_set_number( data, b'heat_trace_power', 500 )
    ssc.data_set_number( data, b'preheat_flux', 50 )
    ssc.data_set_number( data, b'startup_ramp_time', 0 )
    ssc.data_set_number( data, b'startup_target_Tdiff', -5 )
    ssc.data_set_number( data, b'is_rec_startup_from_T_soln', 0 )
    ssc.data_set_number( data, b'is_rec_enforce_min_startup', 0 )
    ssc.data_set_number( data, b'csp.pt.tes.init_hot_htf_percent', 30 )
    ssc.data_set_number( data, b'h_tank', 12 )
    ssc.data_set_number( data, b'cold_tank_max_heat', 15 )
    ssc.data_set_number( data, b'u_tank', 0.40000000000000002 )
    ssc.data_set_number( data, b'tank_pairs', 1 )
    ssc.data_set_number( data, b'cold_tank_Thtr', 280 )
    ssc.data_set_number( data, b'h_tank_min', 1 )
    ssc.data_set_number( data, b'hot_tank_Thtr', 500 )
    ssc.data_set_number( data, b'hot_tank_max_heat', 30 )
    ssc.data_set_number( data, b'tanks_in_parallel', 1 )
    ssc.data_set_number( data, b'h_ctes_tank_min', 1 )
    ssc.data_set_number( data, b'ctes_tshours', 15 )
    ssc.data_set_number( data, b'ctes_field_fl', 4 )
    ssc.data_set_number( data, b'h_ctes_tank', 30 )
    ssc.data_set_number( data, b'u_ctes_tank', 0.40000000000000002 )
    ssc.data_set_number( data, b'ctes_tankpairs', 1 )
    ssc.data_set_number( data, b'T_ctes_cold_design', 5 )
    ssc.data_set_number( data, b'T_ctes_warm_design', 10 )
    ssc.data_set_number( data, b'T_ctes_warm_ini', 20 )
    ssc.data_set_number( data, b'T_ctes_cold_ini', 10 )
    ssc.data_set_number( data, b'f_ctes_warm_ini', 0 )
    ssc.data_set_number( data, b'rad_multiplier', 1.5 )
    ssc.data_set_number( data, b'm_dot_radpanel', 8 )
    ssc.data_set_number( data, b'n_rad_tubes', 100 )
    ssc.data_set_number( data, b'W_rad_tubes', 0.050000000000000003 )
    ssc.data_set_number( data, b'L_rad', 100 )
    ssc.data_set_number( data, b'th_rad_panel', 0.002 )
    ssc.data_set_number( data, b'D_rad_tubes', 0.02 )
    ssc.data_set_number( data, b'k_panel', 235 )
    ssc.data_set_number( data, b'epsilon_radtop', 0.94999999999999996 )
    ssc.data_set_number( data, b'epsilon_radbot', 0.070000000000000007 )
    ssc.data_set_number( data, b'epsilon_radgrnd', 0.90000000000000002 )
    ssc.data_set_number( data, b'L_rad_sections', 10 )
    ssc.data_set_number( data, b'epsilon_radHX', 0.80000000000000004 )
    ssc.data_set_number( data, b'ctes_type', 0 )
    ssc.data_set_number( data, b'helio_area_tot', 1269054.5 )
    ssc.data_set_number( data, b'radiator_unitcost', 13 )
    ssc.data_set_number( data, b'radiator_installcost', 22 )
    ssc.data_set_number( data, b'radiator_fluidcost', 0.34000000357627869 )
    ssc.data_set_number( data, b'radfluid_vol_ratio', 3 )
    ssc.data_set_number( data, b'ctes_cost', 0.69999998807907104 )
    ssc.data_set_number( data, b'rad_pressuredrop', 75 )
    ssc.data_set_number( data, b'pc_config', 1 )
    ssc.data_set_number( data, b'pb_pump_coef', 0.55000000000000004 )
    ssc.data_set_number( data, b'startup_time', 0.5 )
    ssc.data_set_number( data, b'startup_frac', 0.5 )
    ssc.data_set_number( data, b'cycle_max_frac', 1.05 )
    ssc.data_set_number( data, b'cycle_cutoff_frac', 0.20000000000000001 )
    ssc.data_set_number( data, b'q_sby_frac', 0.20000000000000001 )
    ssc.data_set_number( data, b'dT_cw_ref', 10 )
    ssc.data_set_number( data, b'T_amb_des', 42 )
    ssc.data_set_number( data, b'P_boil', 100 )
    ssc.data_set_number( data, b'CT', 2 )
    ssc.data_set_number( data, b'T_approach', 5 )
    ssc.data_set_number( data, b'T_ITD_des', 16 )
    ssc.data_set_number( data, b'P_cond_ratio', 1.0027999999999999 )
    ssc.data_set_number( data, b'pb_bd_frac', 0.02 )
    ssc.data_set_number( data, b'P_cond_min', 2 )
    ssc.data_set_number( data, b'n_pl_inc', 8 )
    F_wc =[ 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ssc.data_set_array( data, b'F_wc',  F_wc)
    ssc.data_set_number( data, b'tech_type', 1 )
    ssc.data_set_number( data, b'ud_f_W_dot_cool_des', 0 )
    ssc.data_set_number( data, b'ud_m_dot_water_cool_des', 0 )
    ssc.data_set_matrix_from_csv( data, b'ud_ind_od', ud_ind_od)
    ssc.data_set_number( data, b'pb_fixed_par', 0.0054999999999999997 )
    ssc.data_set_number( data, b'aux_par', 0.023 )
    ssc.data_set_number( data, b'aux_par_f', 1 )
    ssc.data_set_number( data, b'aux_par_0', 0.48299999999999998 )
    ssc.data_set_number( data, b'aux_par_1', 0.57099999999999995 )
    ssc.data_set_number( data, b'aux_par_2', 0 )
    ssc.data_set_number( data, b'bop_par', 0 )
    ssc.data_set_number( data, b'bop_par_f', 1 )
    ssc.data_set_number( data, b'bop_par_0', 0 )
    ssc.data_set_number( data, b'bop_par_1', 0.48299999999999998 )
    ssc.data_set_number( data, b'bop_par_2', 0 )
    f_turb_tou_periods =[ 1.05, 1, 1, 1, 1, 1, 1, 1, 1 ]
    ssc.data_set_array( data, b'f_turb_tou_periods',  f_turb_tou_periods)
    weekday_schedule = [[ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ]];
    ssc.data_set_matrix( data, b'weekday_schedule', weekday_schedule )
    weekend_schedule = [[ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                        [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                        [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ]];
    ssc.data_set_matrix( data, b'weekend_schedule', weekend_schedule )
    ssc.data_set_number( data, b'is_tod_pc_target_also_pc_max', 0 )
    ssc.data_set_number( data, b'is_dispatch', 0 )
    ssc.data_set_number( data, b'disp_horizon', 48 )
    ssc.data_set_number( data, b'disp_frequency', 24 )
    ssc.data_set_number( data, b'disp_max_iter', 35000 )
    ssc.data_set_number( data, b'disp_timeout', 5 )
    ssc.data_set_number( data, b'disp_mip_gap', 0.001 )
    ssc.data_set_number( data, b'disp_time_weighting', 0.98999999999999999 )
    ssc.data_set_number( data, b'disp_rsu_cost', 950 )
    ssc.data_set_number( data, b'disp_csu_cost', 10000 )
    ssc.data_set_number( data, b'disp_pen_delta_w', 0.10000000000000001 )
    ssc.data_set_number( data, b'disp_inventory_incentive', 0.14999999999999999 )
    ssc.data_set_number( data, b'is_wlim_series', 1 )
    ssc.data_set_array_from_csv( data, b'wlim_series', wlim_series)
    dispatch_sched_weekday = [[ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ]];
    ssc.data_set_matrix( data, b'dispatch_sched_weekday', dispatch_sched_weekday )
    dispatch_sched_weekend = [[ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                              [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], 
                              [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ]];
    ssc.data_set_matrix( data, b'dispatch_sched_weekend', dispatch_sched_weekend )
    ssc.data_set_number( data, b'dispatch_factor1', 2.0640000000000001 )
    ssc.data_set_number( data, b'dispatch_factor2', 1.2 )
    ssc.data_set_number( data, b'dispatch_factor3', 1 )
    ssc.data_set_number( data, b'dispatch_factor4', 1.1000000000000001 )
    ssc.data_set_number( data, b'dispatch_factor5', 0.80000000000000004 )
    ssc.data_set_number( data, b'dispatch_factor6', 0.69999999999999996 )
    ssc.data_set_number( data, b'dispatch_factor7', 1 )
    ssc.data_set_number( data, b'dispatch_factor8', 1 )
    ssc.data_set_number( data, b'dispatch_factor9', 1 )
    ssc.data_set_number( data, b'is_dispatch_series', 0 )
    dispatch_series =[ 0 ]
    ssc.data_set_array( data, b'dispatch_series',  dispatch_series)
    ssc.data_set_number( data, b'rec_height', 21.602900000000002 )
    ssc.data_set_number( data, b'D_rec', 17.649999999999999 )
    ssc.data_set_number( data, b'h_tower', 193.458 )
    ssc.data_set_matrix_from_csv( data, b'helio_positions', helio_positions)
    ssc.data_set_number( data, b'land_area_base', 1847.0377197265625 )
    ssc.data_set_number( data, b'const_per_interest_rate1', 4 )
    ssc.data_set_number( data, b'const_per_interest_rate2', 0 )
    ssc.data_set_number( data, b'const_per_interest_rate3', 0 )
    ssc.data_set_number( data, b'const_per_interest_rate4', 0 )
    ssc.data_set_number( data, b'const_per_interest_rate5', 0 )
    ssc.data_set_number( data, b'const_per_months1', 24 )
    ssc.data_set_number( data, b'const_per_months2', 0 )
    ssc.data_set_number( data, b'const_per_months3', 0 )
    ssc.data_set_number( data, b'const_per_months4', 0 )
    ssc.data_set_number( data, b'const_per_months5', 0 )
    ssc.data_set_number( data, b'const_per_percent1', 100 )
    ssc.data_set_number( data, b'const_per_percent2', 0 )
    ssc.data_set_number( data, b'const_per_percent3', 0 )
    ssc.data_set_number( data, b'const_per_percent4', 0 )
    ssc.data_set_number( data, b'const_per_percent5', 0 )
    ssc.data_set_number( data, b'const_per_upfront_rate1', 1 )
    ssc.data_set_number( data, b'const_per_upfront_rate2', 0 )
    ssc.data_set_number( data, b'const_per_upfront_rate3', 0 )
    ssc.data_set_number( data, b'const_per_upfront_rate4', 0 )
    ssc.data_set_number( data, b'const_per_upfront_rate5', 0 )
    ssc.data_set_number( data, b'adjust:constant', 4 )
    ssc.data_set_number( data, b'sf_adjust:constant', 0 )
    module = ssc.module_create(module_name)    
    ssc.module_exec_set_print( 0 )
    if ssc.module_exec(module, data) == 0:
        print ('nuclear_tes simulation error')
        idx = 1
        msg = ssc.module_log(module, 0)
        while (msg != None):
            print ('    : ' + msg.decode("utf - 8"))
            msg = ssc.module_log(module, idx)
            idx = idx + 1
        SystemExit( "Simulation Error" )
    ssc.module_free(module)
    ssc.data_set_number( data, b'system_use_lifetime_output', 0 )
    ssc.data_set_number( data, b'analysis_period', 25 )
    ssc.data_set_number( data, b'enable_interconnection_limit', 0 )
    ssc.data_set_number( data, b'grid_interconnection_limit_kwac', 20000 )
    ssc.data_set_array_from_csv( data, b'grid_curtailment', grid_curtailment)
    module = ssc.module_create(b'grid')    
    ssc.module_exec_set_print( 0 )
    if ssc.module_exec(module, data) == 0:
        print ('grid simulation error')
        idx = 1
        msg = ssc.module_log(module, 0)
        while (msg != None):
            print ('    : ' + msg.decode("utf - 8"))
            msg = ssc.module_log(module, idx)
            idx = idx + 1
        SystemExit( "Simulation Error" )
    ssc.module_free(module)
    ssc.data_set_number( data, b'ppa_soln_mode', 0 )
    ppa_price_input =[ 0.13 ]
    ssc.data_set_array( data, b'ppa_price_input',  ppa_price_input)
    ssc.data_set_number( data, b'ppa_escalation', 1 )
    federal_tax_rate =[ 21 ]
    ssc.data_set_array( data, b'federal_tax_rate',  federal_tax_rate)
    state_tax_rate =[ 7 ]
    ssc.data_set_array( data, b'state_tax_rate',  state_tax_rate)
    ssc.data_set_number( data, b'property_tax_rate', 0 )
    ssc.data_set_number( data, b'prop_tax_cost_assessed_percent', 100 )
    ssc.data_set_number( data, b'prop_tax_assessed_decline', 0 )
    ssc.data_set_number( data, b'real_discount_rate', 6.4000000000000004 )
    ssc.data_set_number( data, b'inflation_rate', 2.5 )
    ssc.data_set_number( data, b'insurance_rate', 0.5 )
    ssc.data_set_number( data, b'system_capacity', 103500 )
    om_fixed =[ 0 ]
    ssc.data_set_array( data, b'om_fixed',  om_fixed)
    ssc.data_set_number( data, b'om_fixed_escal', 0 )
    om_production =[ 3.5 ]
    ssc.data_set_array( data, b'om_production',  om_production)
    ssc.data_set_number( data, b'om_production_escal', 0 )
    om_capacity =[ 66 ]
    ssc.data_set_array( data, b'om_capacity',  om_capacity)
    ssc.data_set_number( data, b'om_capacity_escal', 0 )
    om_fuel_cost =[ 0 ]
    ssc.data_set_array( data, b'om_fuel_cost',  om_fuel_cost)
    ssc.data_set_number( data, b'om_fuel_cost_escal', 0 )
    om_replacement_cost1 =[ 0 ]
    ssc.data_set_array( data, b'om_replacement_cost1',  om_replacement_cost1)
    ssc.data_set_number( data, b'om_replacement_cost_escal', 0 )
    ssc.data_set_number( data, b'reserves_interest', 1.75 )
    ssc.data_set_number( data, b'equip1_reserve_cost', 0 )
    ssc.data_set_number( data, b'equip1_reserve_freq', 12 )
    ssc.data_set_number( data, b'equip2_reserve_cost', 0 )
    ssc.data_set_number( data, b'equip2_reserve_freq', 15 )
    ssc.data_set_number( data, b'equip3_reserve_cost', 0 )
    ssc.data_set_number( data, b'equip3_reserve_freq', 3 )
    ssc.data_set_number( data, b'equip_reserve_depr_sta', 0 )
    ssc.data_set_number( data, b'equip_reserve_depr_fed', 0 )
    ssc.data_set_number( data, b'itc_fed_amount', 0 )
    ssc.data_set_number( data, b'itc_fed_amount_deprbas_fed', 1 )
    ssc.data_set_number( data, b'itc_fed_amount_deprbas_sta', 1 )
    ssc.data_set_number( data, b'itc_sta_amount', 0 )
    ssc.data_set_number( data, b'itc_sta_amount_deprbas_fed', 0 )
    ssc.data_set_number( data, b'itc_sta_amount_deprbas_sta', 0 )
    ssc.data_set_number( data, b'itc_fed_percent', 26 )
    ssc.data_set_number( data, b'itc_fed_percent_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'itc_fed_percent_deprbas_fed', 1 )
    ssc.data_set_number( data, b'itc_fed_percent_deprbas_sta', 1 )
    ssc.data_set_number( data, b'itc_sta_percent', 0 )
    ssc.data_set_number( data, b'itc_sta_percent_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'itc_sta_percent_deprbas_fed', 0 )
    ssc.data_set_number( data, b'itc_sta_percent_deprbas_sta', 0 )
    ptc_fed_amount =[ 0 ]
    ssc.data_set_array( data, b'ptc_fed_amount',  ptc_fed_amount)
    ssc.data_set_number( data, b'ptc_fed_term', 10 )
    ssc.data_set_number( data, b'ptc_fed_escal', 0 )
    ptc_sta_amount =[ 0 ]
    ssc.data_set_array( data, b'ptc_sta_amount',  ptc_sta_amount)
    ssc.data_set_number( data, b'ptc_sta_term', 10 )
    ssc.data_set_number( data, b'ptc_sta_escal', 0 )
    ssc.data_set_number( data, b'depr_alloc_macrs_5_percent', 90 )
    ssc.data_set_number( data, b'depr_alloc_macrs_15_percent', 1.5 )
    ssc.data_set_number( data, b'depr_alloc_sl_5_percent', 0 )
    ssc.data_set_number( data, b'depr_alloc_sl_15_percent', 2.5 )
    ssc.data_set_number( data, b'depr_alloc_sl_20_percent', 3 )
    ssc.data_set_number( data, b'depr_alloc_sl_39_percent', 0 )
    ssc.data_set_number( data, b'depr_alloc_custom_percent', 0 )
    depr_custom_schedule =[ 0 ]
    ssc.data_set_array( data, b'depr_custom_schedule',  depr_custom_schedule)
    ssc.data_set_number( data, b'depr_bonus_sta', 0 )
    ssc.data_set_number( data, b'depr_bonus_sta_macrs_5', 1 )
    ssc.data_set_number( data, b'depr_bonus_sta_macrs_15', 1 )
    ssc.data_set_number( data, b'depr_bonus_sta_sl_5', 0 )
    ssc.data_set_number( data, b'depr_bonus_sta_sl_15', 0 )
    ssc.data_set_number( data, b'depr_bonus_sta_sl_20', 0 )
    ssc.data_set_number( data, b'depr_bonus_sta_sl_39', 0 )
    ssc.data_set_number( data, b'depr_bonus_sta_custom', 0 )
    ssc.data_set_number( data, b'depr_bonus_fed', 0 )
    ssc.data_set_number( data, b'depr_bonus_fed_macrs_5', 1 )
    ssc.data_set_number( data, b'depr_bonus_fed_macrs_15', 1 )
    ssc.data_set_number( data, b'depr_bonus_fed_sl_5', 0 )
    ssc.data_set_number( data, b'depr_bonus_fed_sl_15', 0 )
    ssc.data_set_number( data, b'depr_bonus_fed_sl_20', 0 )
    ssc.data_set_number( data, b'depr_bonus_fed_sl_39', 0 )
    ssc.data_set_number( data, b'depr_bonus_fed_custom', 0 )
    ssc.data_set_number( data, b'depr_itc_sta_macrs_5', 1 )
    ssc.data_set_number( data, b'depr_itc_sta_macrs_15', 0 )
    ssc.data_set_number( data, b'depr_itc_sta_sl_5', 0 )
    ssc.data_set_number( data, b'depr_itc_sta_sl_15', 0 )
    ssc.data_set_number( data, b'depr_itc_sta_sl_20', 0 )
    ssc.data_set_number( data, b'depr_itc_sta_sl_39', 0 )
    ssc.data_set_number( data, b'depr_itc_sta_custom', 0 )
    ssc.data_set_number( data, b'depr_itc_fed_macrs_5', 1 )
    ssc.data_set_number( data, b'depr_itc_fed_macrs_15', 0 )
    ssc.data_set_number( data, b'depr_itc_fed_sl_5', 0 )
    ssc.data_set_number( data, b'depr_itc_fed_sl_15', 0 )
    ssc.data_set_number( data, b'depr_itc_fed_sl_20', 0 )
    ssc.data_set_number( data, b'depr_itc_fed_sl_39', 0 )
    ssc.data_set_number( data, b'depr_itc_fed_custom', 0 )
    ssc.data_set_number( data, b'ibi_fed_amount', 0 )
    ssc.data_set_number( data, b'ibi_fed_amount_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_fed_amount_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_fed_amount_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_fed_amount_deprbas_sta', 0 )
    ssc.data_set_number( data, b'ibi_sta_amount', 0 )
    ssc.data_set_number( data, b'ibi_sta_amount_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_sta_amount_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_sta_amount_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_sta_amount_deprbas_sta', 0 )
    ssc.data_set_number( data, b'ibi_uti_amount', 0 )
    ssc.data_set_number( data, b'ibi_uti_amount_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_uti_amount_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_uti_amount_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_uti_amount_deprbas_sta', 0 )
    ssc.data_set_number( data, b'ibi_oth_amount', 0 )
    ssc.data_set_number( data, b'ibi_oth_amount_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_oth_amount_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_oth_amount_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_oth_amount_deprbas_sta', 0 )
    ssc.data_set_number( data, b'ibi_fed_percent', 0 )
    ssc.data_set_number( data, b'ibi_fed_percent_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'ibi_fed_percent_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_fed_percent_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_fed_percent_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_fed_percent_deprbas_sta', 0 )
    ssc.data_set_number( data, b'ibi_sta_percent', 0 )
    ssc.data_set_number( data, b'ibi_sta_percent_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'ibi_sta_percent_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_sta_percent_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_sta_percent_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_sta_percent_deprbas_sta', 0 )
    ssc.data_set_number( data, b'ibi_uti_percent', 0 )
    ssc.data_set_number( data, b'ibi_uti_percent_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'ibi_uti_percent_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_uti_percent_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_uti_percent_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_uti_percent_deprbas_sta', 0 )
    ssc.data_set_number( data, b'ibi_oth_percent', 0 )
    ssc.data_set_number( data, b'ibi_oth_percent_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'ibi_oth_percent_tax_fed', 1 )
    ssc.data_set_number( data, b'ibi_oth_percent_tax_sta', 1 )
    ssc.data_set_number( data, b'ibi_oth_percent_deprbas_fed', 0 )
    ssc.data_set_number( data, b'ibi_oth_percent_deprbas_sta', 0 )
    ssc.data_set_number( data, b'cbi_fed_amount', 0 )
    ssc.data_set_number( data, b'cbi_fed_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'cbi_fed_tax_fed', 1 )
    ssc.data_set_number( data, b'cbi_fed_tax_sta', 1 )
    ssc.data_set_number( data, b'cbi_fed_deprbas_fed', 0 )
    ssc.data_set_number( data, b'cbi_fed_deprbas_sta', 0 )
    ssc.data_set_number( data, b'cbi_sta_amount', 0 )
    ssc.data_set_number( data, b'cbi_sta_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'cbi_sta_tax_fed', 1 )
    ssc.data_set_number( data, b'cbi_sta_tax_sta', 1 )
    ssc.data_set_number( data, b'cbi_sta_deprbas_fed', 0 )
    ssc.data_set_number( data, b'cbi_sta_deprbas_sta', 0 )
    ssc.data_set_number( data, b'cbi_uti_amount', 0 )
    ssc.data_set_number( data, b'cbi_uti_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'cbi_uti_tax_fed', 1 )
    ssc.data_set_number( data, b'cbi_uti_tax_sta', 1 )
    ssc.data_set_number( data, b'cbi_uti_deprbas_fed', 0 )
    ssc.data_set_number( data, b'cbi_uti_deprbas_sta', 0 )
    ssc.data_set_number( data, b'cbi_oth_amount', 0 )
    ssc.data_set_number( data, b'cbi_oth_maxvalue', 9.9999999999999998e+37 )
    ssc.data_set_number( data, b'cbi_oth_tax_fed', 1 )
    ssc.data_set_number( data, b'cbi_oth_tax_sta', 1 )
    ssc.data_set_number( data, b'cbi_oth_deprbas_fed', 0 )
    ssc.data_set_number( data, b'cbi_oth_deprbas_sta', 0 )
    pbi_fed_amount =[ 0 ]
    ssc.data_set_array( data, b'pbi_fed_amount',  pbi_fed_amount)
    ssc.data_set_number( data, b'pbi_fed_term', 0 )
    ssc.data_set_number( data, b'pbi_fed_escal', 0 )
    ssc.data_set_number( data, b'pbi_fed_tax_fed', 1 )
    ssc.data_set_number( data, b'pbi_fed_tax_sta', 1 )
    pbi_sta_amount =[ 0 ]
    ssc.data_set_array( data, b'pbi_sta_amount',  pbi_sta_amount)
    ssc.data_set_number( data, b'pbi_sta_term', 0 )
    ssc.data_set_number( data, b'pbi_sta_escal', 0 )
    ssc.data_set_number( data, b'pbi_sta_tax_fed', 1 )
    ssc.data_set_number( data, b'pbi_sta_tax_sta', 1 )
    pbi_uti_amount =[ 0 ]
    ssc.data_set_array( data, b'pbi_uti_amount',  pbi_uti_amount)
    ssc.data_set_number( data, b'pbi_uti_term', 0 )
    ssc.data_set_number( data, b'pbi_uti_escal', 0 )
    ssc.data_set_number( data, b'pbi_uti_tax_fed', 1 )
    ssc.data_set_number( data, b'pbi_uti_tax_sta', 1 )
    pbi_oth_amount =[ 0 ]
    ssc.data_set_array( data, b'pbi_oth_amount',  pbi_oth_amount)
    ssc.data_set_number( data, b'pbi_oth_term', 0 )
    ssc.data_set_number( data, b'pbi_oth_escal', 0 )
    ssc.data_set_number( data, b'pbi_oth_tax_fed', 1 )
    ssc.data_set_number( data, b'pbi_oth_tax_sta', 1 )
    ssc.data_set_number( data, b'term_tenor', 18 )
    ssc.data_set_number( data, b'term_int_rate', 7 )
    ssc.data_set_number( data, b'dscr', 1.3 )
    ssc.data_set_number( data, b'dscr_reserve_months', 6 )
    ssc.data_set_number( data, b'debt_percent', 50 )
    ssc.data_set_number( data, b'debt_option', 1 )
    ssc.data_set_number( data, b'payment_option', 0 )
    ssc.data_set_number( data, b'cost_debt_closing', 450000 )
    ssc.data_set_number( data, b'cost_debt_fee', 2.75 )
    ssc.data_set_number( data, b'months_working_reserve', 6 )
    ssc.data_set_number( data, b'months_receivables_reserve', 0 )
    ssc.data_set_number( data, b'cost_other_financing', 0 )
    ssc.data_set_number( data, b'flip_target_percent', 11 )
    ssc.data_set_number( data, b'flip_target_year', 20 )
    ssc.data_set_number( data, b'pbi_fed_for_ds', 0 )
    ssc.data_set_number( data, b'pbi_sta_for_ds', 0 )
    ssc.data_set_number( data, b'pbi_uti_for_ds', 0 )
    ssc.data_set_number( data, b'pbi_oth_for_ds', 0 )
    degradation =[ 0 ]
    ssc.data_set_array( data, b'degradation',  degradation)
    roe_input =[ 0 ]
    ssc.data_set_array( data, b'roe_input',  roe_input)
    ssc.data_set_number( data, b'loan_moratorium', 0 )
    ssc.data_set_number( data, b'system_use_recapitalization', 0 )
    ssc.data_set_number( data, b'total_installed_cost', 673465472 )
    ssc.data_set_number( data, b'salvage_percentage', 0 )
    ssc.data_set_number( data, b'construction_financing_cost', 33673272 )
    ssc.data_set_number( data, b'depr_stabas_method', 1 )
    ssc.data_set_number( data, b'depr_fedbas_method', 1 )
    ssc.data_set_number( data, b'cp_capacity_payment_esc', 0 )
    ssc.data_set_number( data, b'cp_capacity_payment_type', 0 )
    cp_capacity_payment_amount =[ 0 ]
    ssc.data_set_array( data, b'cp_capacity_payment_amount',  cp_capacity_payment_amount)
    cp_capacity_credit_percent =[ 0 ]
    ssc.data_set_array( data, b'cp_capacity_credit_percent',  cp_capacity_credit_percent)
    ssc.data_set_number( data, b'cp_system_nameplate', 103.5 )
    ssc.data_set_number( data, b'cp_battery_nameplate', 0 )
    grid_curtailment_price =[ 0 ]
    ssc.data_set_array( data, b'grid_curtailment_price',  grid_curtailment_price)
    ssc.data_set_number( data, b'grid_curtailment_price_esc', 0 )
    module = ssc.module_create(b'singleowner')    
    ssc.module_exec_set_print( 0 )
    if ssc.module_exec(module, data) == 0:
        print ('singleowner simulation error')
        idx = 1
        msg = ssc.module_log(module, 0)
        while (msg != None):
            print ('    : ' + msg.decode("utf - 8"))
            msg = ssc.module_log(module, idx)
            idx = idx + 1
        SystemExit( "Simulation Error" )
    ssc.module_free(module)
    print('                        CSP')
    annual_energy = ssc.data_get_number(data, b'annual_energy') / 1e6
    print ('Annual energy (year 1)        = ', annual_energy, ' GWh')
    capacity_factor = ssc.data_get_number(data, b'capacity_factor')
    print ('Capacity factor (year 1)      = ', capacity_factor, ' %')
    annual_total_water_use = ssc.data_get_number(data, b'annual_total_water_use')
    print ('Annual Water Usage            = ', annual_total_water_use, ' m3')
    ppa = ssc.data_get_number(data, b'ppa')
    print ('PPA price (year 1)            = ', ppa)
    lppa_nom = ssc.data_get_number(data, b'lppa_nom')
    print ('Levelized PPA price (nominal) = ', lppa_nom)
    lppa_real = ssc.data_get_number(data, b'lppa_real')
    print ('Levelized PPA price (real)    = ', lppa_real)
    lcoe_nom = ssc.data_get_number(data, b'lcoe_nom')
    print ('Levelized COE (nominal)       = ', lcoe_nom)
    lcoe_real = ssc.data_get_number(data, b'lcoe_real')
    print ('Levelized COE (real)          = ', lcoe_real)
    project_return_aftertax_npv = ssc.data_get_number(data, b'project_return_aftertax_npv')
    print ('Net present value             = ', project_return_aftertax_npv)
    flip_actual_irr = ssc.data_get_number(data, b'flip_actual_irr')
    print ('Internal rate of return (IRR) = ', flip_actual_irr)
    flip_actual_year = ssc.data_get_number(data, b'flip_actual_year')
    print ('Year IRR is achieved          = ', flip_actual_year)
    project_return_aftertax_irr = ssc.data_get_number(data, b'project_return_aftertax_irr')
    print ('IRR at end of project         = ', project_return_aftertax_irr)
    cost_installed = ssc.data_get_number(data, b'cost_installed')
    print ('Net capital cost              = ', cost_installed)
    size_of_equity = ssc.data_get_number(data, b'size_of_equity')
    print ('Equity                        = ', size_of_equity)
    size_of_debt = ssc.data_get_number(data, b'size_of_debt')
    print ('Size of debt                  = ', size_of_debt)
    a_sf = ssc.data_get_number(data, b'A_sf')
    print ('Area of solar field           = ', a_sf)


# =============================================================================
#     Plotting
# =============================================================================

    def get_array(out_str):
        out_array = ssc.data_get_array(data,out_str.encode('utf-8'))
        out_array = np.asarray(out_array)[0:hours]
        return out_array
        
    p_cycle       = get_array('P_cycle')
    gen           = get_array('gen') / 1e3
    p_cool        = get_array('P_cooling_tower_tot')
    q_dot_rec_in  = get_array('q_dot_rec_inc')
    m_dot         = get_array('m_dot_rec')
    T_pc_in       = get_array('T_pc_in')
    T_pc_out      = get_array('T_pc_out')
    e_ch_tes      = get_array('e_ch_tes')
    t_plot        = get_array('time_hr') / 24
    op_mode_1      = get_array('op_mode_1')
    
    #operating modes, copied from ssc/tcs/csp_solver.cpp
    n_modes, modes_order = np.unique(op_mode_1,return_index=True)
    n_modes = n_modes[np.argsort(modes_order)] # re-order modes by first appearance of each
    op_modes_list = [
        "ENTRY_MODE",
        "CR_OFF__PC_OFF__TES_OFF,"
        "CR_SU__PC_OFF__TES_OFF",
        "CR_ON__PC_SU__TES_OFF",
        "CR_ON__PC_SB__TES_OFF",        
        "CR_ON__PC_RM_HI__TES_OFF",
        "CR_ON__PC_RM_LO__TES_OFF",        
        "CR_DF__PC_MAX__TES_OFF",
        "CR_OFF__PC_SU__TES_DC",
        "CR_ON__PC_OFF__TES_CH",
        "SKIP_10",
        "CR_ON__PC_TARGET__TES_CH",
        "CR_ON__PC_TARGET__TES_DC",
        "CR_ON__PC_RM_LO__TES_EMPTY",
        "CR_DF__PC_OFF__TES_FULL",       
        "CR_OFF__PC_SB__TES_DC",
        "CR_OFF__PC_MIN__TES_EMPTY",
        "CR_OFF__PC_RM_LO__TES_EMPTY",
        "CR_ON__PC_SB__TES_CH",
        "CR_SU__PC_MIN__TES_EMPTY",
        "SKIP_20",
        "CR_SU__PC_SB__TES_DC",
        "CR_ON__PC_SB__TES_DC",
        "CR_OFF__PC_TARGET__TES_DC",
        "CR_SU__PC_TARGET__TES_DC",
        "CR_ON__PC_RM_HI__TES_FULL",
        "CR_ON__PC_MIN__TES_EMPTY",
        "CR_SU__PC_RM_LO__TES_EMPTY",
        "CR_DF__PC_MAX__TES_FULL",
        "CR_ON__PC_SB__TES_FULL",
        "SKIP_30",
        "CR_SU__PC_SU__TES_DC",
        "CR_ON__PC_SU__TES_CH",
        "CR_DF__PC_SU__TES_FULL",
        "CR_DF__PC_SU__TES_OFF",
        "CR_TO_COLD__PC_TARGET__TES_DC",
        "CR_TO_COLD__PC_RM_LO__TES_EMPTY",
        "CR_TO_COLD__PC_SB__TES_DC",
        "CR_TO_COLD__PC_MIN__TES_EMPTY",
        "CR_TO_COLD__PC_OFF__TES_OFF",
        "SKIP_40",
        "CR_TO_COLD__PC_SU__TES_DC" ]


    lp = 16 #labelpad
    fs = 12 #fontsize
    lw = 2  #linewidth
    fsl = 'x-small'      #fontsize legend
    loc = 'upper right'  #location of legend
    
    
    fig = plt.figure(figsize=[10,8])
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    # Energy plot
    ax1.plot(t_plot, e_ch_tes, linewidth = lw, label='Salt Charge Level (MWht)')
    ax1.plot(t_plot, p_cycle, linewidth = lw, label='P_cycle (Electric)')
    ax1.plot(t_plot, q_dot_rec_in, linewidth = lw, label='Q_dot to Salt (Thermal)')
    ax1.plot(t_plot, gen, linewidth = lw, label='Power generated')
    ax1.set_ylabel('Power (MW)', labelpad=lp, fontsize=fs, fontweight='bold')
    ax1.legend(loc=loc,fontsize=fsl)
    
    # Mass flow rate plot
    ax2.plot(t_plot, m_dot, linewidth = lw, label='m_dot_water_pc')
    ax2.set_ylabel('Mass flow (kg/s)', labelpad=lp, fontsize=fs, fontweight='bold')
    ax2.legend(loc=loc,fontsize=fsl)
    
    # Temperature plot
    ax3.plot(t_plot, T_pc_in, linewidth = lw, label='PC HTF inlet')
    ax3.plot(t_plot, T_pc_out, linewidth = lw, label='PC HTF outlet')
    ax3.set_xlabel('Time (days)', labelpad=lp, fontsize=fs, fontweight='bold')
    ax3.set_ylabel('Temperature (C)', labelpad=lp, fontsize=fs, fontweight='bold')
    ax3.legend(loc=loc,fontsize=fsl)
    
    # operating modes time history
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(t_plot, op_mode_1)
    for op in n_modes:
        inds = op_mode_1 == op
        ax.plot(t_plot[inds], op_mode_1[inds], '.', label=op_modes_list[int(op)])
    ax.legend(loc='best')
