from PySSC import PySSC

def battrun(weather, battcapacity):
	ssc = PySSC()
# 	print ('Current folder = C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep')
# 	print ('SSC Version = ', ssc.version())
# 	print ('SSC Build Information = ', ssc.build_info().decode("utf - 8"))
	ssc.module_exec_set_print(0)
	data = ssc.data_create()
	ssc.data_set_string( data, b'solar_resource_file', weather );
	ssc.data_set_number( data, b'transformer_no_load_loss', 0 )
	ssc.data_set_number( data, b'transformer_load_loss', 0 )
	ssc.data_set_number( data, b'system_use_lifetime_output', 1 )
	ssc.data_set_number( data, b'analysis_period', 25 )
	dc_degradation =[ 0.5 ];
	ssc.data_set_array( data, b'dc_degradation',  dc_degradation);
	ssc.data_set_number( data, b'en_dc_lifetime_losses', 0 )
	ssc.data_set_array_from_csv( data, b'dc_lifetime_losses', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/dc_lifetime_losses.csv');
	ssc.data_set_number( data, b'en_ac_lifetime_losses', 0 )
	ssc.data_set_array_from_csv( data, b'ac_lifetime_losses', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/ac_lifetime_losses.csv');
	ssc.data_set_number( data, b'save_full_lifetime_variables', 1 )
	ssc.data_set_number( data, b'en_snow_model', 0 )
	ssc.data_set_number( data, b'system_capacity', 498.71957397460938 )
	ssc.data_set_number( data, b'use_wf_albedo', 0 )
	albedo =[ 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001 ];
	ssc.data_set_array( data, b'albedo',  albedo);
	ssc.data_set_number( data, b'irrad_mode', 0 )
	ssc.data_set_number( data, b'sky_model', 2 )
	ssc.data_set_number( data, b'inverter_count', 7 )
	ssc.data_set_number( data, b'enable_mismatch_vmax_calc', 0 )
	ssc.data_set_number( data, b'subarray1_nstrings', 134 )
	ssc.data_set_number( data, b'subarray1_modules_per_string', 12 )
	ssc.data_set_number( data, b'subarray1_mppt_input', 1 )
	ssc.data_set_number( data, b'subarray1_tilt', 20 )
	ssc.data_set_number( data, b'subarray1_tilt_eq_lat', 0 )
	ssc.data_set_number( data, b'subarray1_azimuth', 180 )
	ssc.data_set_number( data, b'subarray1_track_mode', 0 )
	ssc.data_set_number( data, b'subarray1_rotlim', 45 )
	ssc.data_set_number( data, b'subarray1_shade_mode', 0 )
	ssc.data_set_number( data, b'subarray1_gcr', 0.29999999999999999 )
	subarray1_monthly_tilt =[ 40, 40, 40, 20, 20, 20, 20, 20, 20, 40, 40, 40 ];
	ssc.data_set_array( data, b'subarray1_monthly_tilt',  subarray1_monthly_tilt);
	subarray1_soiling =[ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ];
	ssc.data_set_array( data, b'subarray1_soiling',  subarray1_soiling);
	ssc.data_set_number( data, b'subarray1_rear_irradiance_loss', 0 )
	ssc.data_set_number( data, b'subarray1_mismatch_loss', 2 )
	ssc.data_set_number( data, b'subarray1_diodeconn_loss', 0.5 )
	ssc.data_set_number( data, b'subarray1_dcwiring_loss', 2 )
	ssc.data_set_number( data, b'subarray1_tracking_loss', 0 )
	ssc.data_set_number( data, b'subarray1_nameplate_loss', 0 )
	ssc.data_set_number( data, b'subarray2_rear_irradiance_loss', 0 )
	ssc.data_set_number( data, b'subarray2_mismatch_loss', 2 )
	ssc.data_set_number( data, b'subarray2_diodeconn_loss', 0.5 )
	ssc.data_set_number( data, b'subarray2_dcwiring_loss', 2 )
	ssc.data_set_number( data, b'subarray2_tracking_loss', 0 )
	ssc.data_set_number( data, b'subarray2_nameplate_loss', 0 )
	ssc.data_set_number( data, b'subarray3_rear_irradiance_loss', 0 )
	ssc.data_set_number( data, b'subarray3_mismatch_loss', 2 )
	ssc.data_set_number( data, b'subarray3_diodeconn_loss', 0.5 )
	ssc.data_set_number( data, b'subarray3_dcwiring_loss', 2 )
	ssc.data_set_number( data, b'subarray3_tracking_loss', 0 )
	ssc.data_set_number( data, b'subarray3_nameplate_loss', 0 )
	ssc.data_set_number( data, b'subarray4_rear_irradiance_loss', 0 )
	ssc.data_set_number( data, b'subarray4_mismatch_loss', 2 )
	ssc.data_set_number( data, b'subarray4_diodeconn_loss', 0.5 )
	ssc.data_set_number( data, b'subarray4_dcwiring_loss', 2 )
	ssc.data_set_number( data, b'subarray4_tracking_loss', 0 )
	ssc.data_set_number( data, b'subarray4_nameplate_loss', 0 )
	ssc.data_set_number( data, b'dcoptimizer_loss', 0 )
	ssc.data_set_number( data, b'acwiring_loss', 1 )
	ssc.data_set_number( data, b'transmission_loss', 0 )
	ssc.data_set_number( data, b'subarray1_mod_orient', 0 )
	ssc.data_set_number( data, b'subarray1_nmodx', 7 )
	ssc.data_set_number( data, b'subarray1_nmody', 2 )
	ssc.data_set_number( data, b'subarray1_backtrack', 0 )
	ssc.data_set_number( data, b'subarray2_enable', 0 )
	ssc.data_set_number( data, b'subarray2_modules_per_string', 0 )
	ssc.data_set_number( data, b'subarray2_nstrings', 0 )
	ssc.data_set_number( data, b'subarray2_mppt_input', 1 )
	ssc.data_set_number( data, b'subarray2_tilt', 20 )
	ssc.data_set_number( data, b'subarray2_tilt_eq_lat', 0 )
	ssc.data_set_number( data, b'subarray2_azimuth', 180 )
	ssc.data_set_number( data, b'subarray2_track_mode', 0 )
	ssc.data_set_number( data, b'subarray2_rotlim', 45 )
	ssc.data_set_number( data, b'subarray2_shade_mode', 0 )
	ssc.data_set_number( data, b'subarray2_gcr', 0.29999999999999999 )
	subarray2_monthly_tilt =[ 40, 40, 40, 20, 20, 20, 20, 20, 20, 40, 40, 40 ];
	ssc.data_set_array( data, b'subarray2_monthly_tilt',  subarray2_monthly_tilt);
	subarray2_soiling =[ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ];
	ssc.data_set_array( data, b'subarray2_soiling',  subarray2_soiling);
	ssc.data_set_number( data, b'subarray2_mod_orient', 0 )
	ssc.data_set_number( data, b'subarray2_nmodx', 9 )
	ssc.data_set_number( data, b'subarray2_nmody', 2 )
	ssc.data_set_number( data, b'subarray2_backtrack', 0 )
	ssc.data_set_number( data, b'subarray3_enable', 0 )
	ssc.data_set_number( data, b'subarray3_modules_per_string', 0 )
	ssc.data_set_number( data, b'subarray3_nstrings', 0 )
	ssc.data_set_number( data, b'subarray3_mppt_input', 1 )
	ssc.data_set_number( data, b'subarray3_tilt', 20 )
	ssc.data_set_number( data, b'subarray3_tilt_eq_lat', 0 )
	ssc.data_set_number( data, b'subarray3_azimuth', 180 )
	ssc.data_set_number( data, b'subarray3_track_mode', 0 )
	ssc.data_set_number( data, b'subarray3_rotlim', 45 )
	ssc.data_set_number( data, b'subarray3_shade_mode', 0 )
	ssc.data_set_number( data, b'subarray3_gcr', 0.29999999999999999 )
	subarray3_monthly_tilt =[ 40, 40, 40, 20, 20, 20, 20, 20, 20, 40, 40, 40 ];
	ssc.data_set_array( data, b'subarray3_monthly_tilt',  subarray3_monthly_tilt);
	subarray3_soiling =[ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ];
	ssc.data_set_array( data, b'subarray3_soiling',  subarray3_soiling);
	ssc.data_set_number( data, b'subarray3_mod_orient', 0 )
	ssc.data_set_number( data, b'subarray3_nmodx', 9 )
	ssc.data_set_number( data, b'subarray3_nmody', 2 )
	ssc.data_set_number( data, b'subarray3_backtrack', 0 )
	ssc.data_set_number( data, b'subarray4_enable', 0 )
	ssc.data_set_number( data, b'subarray4_modules_per_string', 0 )
	ssc.data_set_number( data, b'subarray4_nstrings', 0 )
	ssc.data_set_number( data, b'subarray4_mppt_input', 1 )
	ssc.data_set_number( data, b'subarray4_tilt', 20 )
	ssc.data_set_number( data, b'subarray4_tilt_eq_lat', 0 )
	ssc.data_set_number( data, b'subarray4_azimuth', 180 )
	ssc.data_set_number( data, b'subarray4_track_mode', 0 )
	ssc.data_set_number( data, b'subarray4_rotlim', 45 )
	ssc.data_set_number( data, b'subarray4_shade_mode', 0 )
	ssc.data_set_number( data, b'subarray4_gcr', 0.29999999999999999 )
	subarray4_monthly_tilt =[ 40, 40, 40, 20, 20, 20, 20, 20, 20, 40, 40, 40 ];
	ssc.data_set_array( data, b'subarray4_monthly_tilt',  subarray4_monthly_tilt);
	subarray4_soiling =[ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ];
	ssc.data_set_array( data, b'subarray4_soiling',  subarray4_soiling);
	ssc.data_set_number( data, b'subarray4_mod_orient', 0 )
	ssc.data_set_number( data, b'subarray4_nmodx', 9 )
	ssc.data_set_number( data, b'subarray4_nmody', 2 )
	ssc.data_set_number( data, b'subarray4_backtrack', 0 )
	ssc.data_set_number( data, b'module_model', 1 )
	ssc.data_set_number( data, b'module_aspect_ratio', 1.7 )
	ssc.data_set_number( data, b'spe_area', 0.74073999999999995 )
	ssc.data_set_number( data, b'spe_rad0', 200 )
	ssc.data_set_number( data, b'spe_rad1', 400 )
	ssc.data_set_number( data, b'spe_rad2', 600 )
	ssc.data_set_number( data, b'spe_rad3', 800 )
	ssc.data_set_number( data, b'spe_rad4', 1000 )
	ssc.data_set_number( data, b'spe_eff0', 13.5 )
	ssc.data_set_number( data, b'spe_eff1', 13.5 )
	ssc.data_set_number( data, b'spe_eff2', 13.5 )
	ssc.data_set_number( data, b'spe_eff3', 13.5 )
	ssc.data_set_number( data, b'spe_eff4', 13.5 )
	ssc.data_set_number( data, b'spe_reference', 4 )
	ssc.data_set_number( data, b'spe_module_structure', 0 )
	ssc.data_set_number( data, b'spe_a', -3.5600000000000001 )
	ssc.data_set_number( data, b'spe_b', -0.074999999999999997 )
	ssc.data_set_number( data, b'spe_dT', 3 )
	ssc.data_set_number( data, b'spe_temp_coeff', -0.5 )
	ssc.data_set_number( data, b'spe_fd', 1 )
	ssc.data_set_number( data, b'spe_vmp', 30 )
	ssc.data_set_number( data, b'spe_voc', 36 )
	ssc.data_set_number( data, b'spe_is_bifacial', 0 )
	ssc.data_set_number( data, b'spe_bifacial_transmission_factor', 0.012999999999999999 )
	ssc.data_set_number( data, b'spe_bifaciality', 0.65000000000000002 )
	ssc.data_set_number( data, b'spe_bifacial_ground_clearance_height', 1 )
	ssc.data_set_number( data, b'cec_area', 1.631 )
	ssc.data_set_number( data, b'cec_a_ref', 2.5776400000000002 )
	ssc.data_set_number( data, b'cec_adjust', 22.909172999999999 )
	ssc.data_set_number( data, b'cec_alpha_sc', 0.003735 )
	ssc.data_set_number( data, b'cec_beta_oc', -0.175619 )
	ssc.data_set_number( data, b'cec_gamma_r', -0.38600000000000001 )
	ssc.data_set_number( data, b'cec_i_l_ref', 6.0537280000000004 )
	ssc.data_set_number( data, b'cec_i_mp_ref', 5.6699999999999999 )
	ssc.data_set_number( data, b'cec_i_o_ref', 8.3604300000000002e-11 )
	ssc.data_set_number( data, b'cec_i_sc_ref', 6.0499999999999998 )
	ssc.data_set_number( data, b'cec_n_s', 96 )
	ssc.data_set_number( data, b'cec_r_s', 0.30812 )
	ssc.data_set_number( data, b'cec_r_sh_ref', 500.06892099999999 )
	ssc.data_set_number( data, b'cec_t_noct', 46 )
	ssc.data_set_number( data, b'cec_v_mp_ref', 54.700000000000003 )
	ssc.data_set_number( data, b'cec_v_oc_ref', 64.400000000000006 )
	ssc.data_set_number( data, b'cec_temp_corr_mode', 0 )
	ssc.data_set_number( data, b'cec_is_bifacial', 0 )
	ssc.data_set_number( data, b'cec_bifacial_transmission_factor', 0.012999999999999999 )
	ssc.data_set_number( data, b'cec_bifaciality', 0.65000000000000002 )
	ssc.data_set_number( data, b'cec_bifacial_ground_clearance_height', 1 )
	ssc.data_set_number( data, b'cec_standoff', 6 )
	ssc.data_set_number( data, b'cec_height', 0 )
	ssc.data_set_number( data, b'cec_mounting_config', 0 )
	ssc.data_set_number( data, b'cec_heat_transfer', 0 )
	ssc.data_set_number( data, b'cec_mounting_orientation', 0 )
	ssc.data_set_number( data, b'cec_gap_spacing', 0.050000000000000003 )
	ssc.data_set_number( data, b'cec_module_width', 1 )
	ssc.data_set_number( data, b'cec_module_length', 1.6310000419616699 )
	ssc.data_set_number( data, b'cec_array_rows', 1 )
	ssc.data_set_number( data, b'cec_array_cols', 10 )
	ssc.data_set_number( data, b'cec_backside_temp', 20 )
	ssc.data_set_number( data, b'cec_transient_thermal_model_unit_mass', 11.091900000000001 )
	ssc.data_set_number( data, b'6par_celltech', 1 )
	ssc.data_set_number( data, b'6par_vmp', 30 )
	ssc.data_set_number( data, b'6par_imp', 6 )
	ssc.data_set_number( data, b'6par_voc', 37 )
	ssc.data_set_number( data, b'6par_isc', 7 )
	ssc.data_set_number( data, b'6par_bvoc', -0.10999999940395355 )
	ssc.data_set_number( data, b'6par_aisc', 0.0040000001899898052 )
	ssc.data_set_number( data, b'6par_gpmp', -0.40999999999999998 )
	ssc.data_set_number( data, b'6par_nser', 60 )
	ssc.data_set_number( data, b'6par_area', 1.3 )
	ssc.data_set_number( data, b'6par_tnoct', 46 )
	ssc.data_set_number( data, b'6par_standoff', 6 )
	ssc.data_set_number( data, b'6par_mounting', 0 )
	ssc.data_set_number( data, b'6par_is_bifacial', 0 )
	ssc.data_set_number( data, b'6par_bifacial_transmission_factor', 0.012999999999999999 )
	ssc.data_set_number( data, b'6par_bifaciality', 0.65000000000000002 )
	ssc.data_set_number( data, b'6par_bifacial_ground_clearance_height', 1 )
	ssc.data_set_number( data, b'6par_transient_thermal_model_unit_mass', 11.091900000000001 )
	ssc.data_set_number( data, b'snl_module_structure', 0 )
	ssc.data_set_number( data, b'snl_a', -3.6200000000000001 )
	ssc.data_set_number( data, b'snl_b', -0.074999999999999997 )
	ssc.data_set_number( data, b'snl_dtc', 3 )
	ssc.data_set_number( data, b'snl_ref_a', -3.619999885559082 )
	ssc.data_set_number( data, b'snl_ref_b', -0.075000002980232239 )
	ssc.data_set_number( data, b'snl_ref_dT', 3 )
	ssc.data_set_number( data, b'snl_fd', 1 )
	ssc.data_set_number( data, b'snl_a0', 0.94045000000000001 )
	ssc.data_set_number( data, b'snl_a1', 0.052641 )
	ssc.data_set_number( data, b'snl_a2', -0.0093897000000000008 )
	ssc.data_set_number( data, b'snl_a3', 0.00072623000000000002 )
	ssc.data_set_number( data, b'snl_a4', -1.9938000000000001e-05 )
	ssc.data_set_number( data, b'snl_aimp', -0.00038000000000000002 )
	ssc.data_set_number( data, b'snl_aisc', 0.00060999999999999997 )
	ssc.data_set_number( data, b'snl_area', 1.244 )
	ssc.data_set_number( data, b'snl_b0', 1 )
	ssc.data_set_number( data, b'snl_b1', -0.0024380000000000001 )
	ssc.data_set_number( data, b'snl_b2', 0.00031030000000000001 )
	ssc.data_set_number( data, b'snl_b3', -1.2459999999999999e-05 )
	ssc.data_set_number( data, b'snl_b4', 2.11e-07 )
	ssc.data_set_number( data, b'snl_b5', -1.3600000000000001e-09 )
	ssc.data_set_number( data, b'snl_bvmpo', -0.13900000000000001 )
	ssc.data_set_number( data, b'snl_bvoco', -0.13600000000000001 )
	ssc.data_set_number( data, b'snl_c0', 1.0039 )
	ssc.data_set_number( data, b'snl_c1', -0.0038999999999999998 )
	ssc.data_set_number( data, b'snl_c2', 0.29106599999999999 )
	ssc.data_set_number( data, b'snl_c3', -4.7354599999999998 )
	ssc.data_set_number( data, b'snl_c4', 0.99419999999999997 )
	ssc.data_set_number( data, b'snl_c5', 0.0057999999999999996 )
	ssc.data_set_number( data, b'snl_c6', 1.0723 )
	ssc.data_set_number( data, b'snl_c7', -0.072300000000000003 )
	ssc.data_set_number( data, b'snl_impo', 5.25 )
	ssc.data_set_number( data, b'snl_isco', 5.75 )
	ssc.data_set_number( data, b'snl_ixo', 5.6500000000000004 )
	ssc.data_set_number( data, b'snl_ixxo', 3.8500000000000001 )
	ssc.data_set_number( data, b'snl_mbvmp', 0 )
	ssc.data_set_number( data, b'snl_mbvoc', 0 )
	ssc.data_set_number( data, b'snl_n', 1.2210000000000001 )
	ssc.data_set_number( data, b'snl_series_cells', 72 )
	ssc.data_set_number( data, b'snl_vmpo', 40 )
	ssc.data_set_number( data, b'snl_voco', 47.700000000000003 )
	ssc.data_set_number( data, b'snl_transient_thermal_model_unit_mass', 11.091900000000001 )
	ssc.data_set_number( data, b'sd11par_nser', 116 )
	ssc.data_set_number( data, b'sd11par_area', 0.71999999999999997 )
	ssc.data_set_number( data, b'sd11par_AMa0', 0.94169999999999998 )
	ssc.data_set_number( data, b'sd11par_AMa1', 0.065159999999999996 )
	ssc.data_set_number( data, b'sd11par_AMa2', -0.020219999999999998 )
	ssc.data_set_number( data, b'sd11par_AMa3', 0.0021900000000000001 )
	ssc.data_set_number( data, b'sd11par_AMa4', -9.1000000000000003e-05 )
	ssc.data_set_number( data, b'sd11par_glass', 0 )
	ssc.data_set_number( data, b'sd11par_tnoct', 44.899999999999999 )
	ssc.data_set_number( data, b'sd11par_standoff', 6 )
	ssc.data_set_number( data, b'sd11par_mounting', 0 )
	ssc.data_set_number( data, b'sd11par_Vmp0', 64.599998474121094 )
	ssc.data_set_number( data, b'sd11par_Imp0', 1.0499999523162842 )
	ssc.data_set_number( data, b'sd11par_Voc0', 87 )
	ssc.data_set_number( data, b'sd11par_Isc0', 1.1799999475479126 )
	ssc.data_set_number( data, b'sd11par_alphaIsc', 0.000472001 )
	ssc.data_set_number( data, b'sd11par_n', 1.4507099999999999 )
	ssc.data_set_number( data, b'sd11par_Il', 1.1895100000000001 )
	ssc.data_set_number( data, b'sd11par_Io', 2.08522e-09 )
	ssc.data_set_number( data, b'sd11par_Egref', 0.73766799999999999 )
	ssc.data_set_number( data, b'sd11par_d1', 13.5504 )
	ssc.data_set_number( data, b'sd11par_d2', -0.0769735 )
	ssc.data_set_number( data, b'sd11par_d3', 0.23732700000000001 )
	ssc.data_set_number( data, b'sd11par_c1', 1930.1500000000001 )
	ssc.data_set_number( data, b'sd11par_c2', 474.63999999999999 )
	ssc.data_set_number( data, b'sd11par_c3', 1.48746 )
	ssc.data_set_number( data, b'inverter_model', 0 )
	ssc.data_set_number( data, b'mppt_low_inverter', 570 )
	ssc.data_set_number( data, b'mppt_hi_inverter', 800 )
	ssc.data_set_number( data, b'inv_num_mppt', 1 )
	ssc.data_set_number( data, b'inv_snl_c0', -2.0614739999999999e-07 )
	ssc.data_set_number( data, b'inv_snl_c1', 2.6999999999999999e-05 )
	ssc.data_set_number( data, b'inv_snl_c2', 0.0026059999999999998 )
	ssc.data_set_number( data, b'inv_snl_c3', 0.00050100000000000003 )
	ssc.data_set_number( data, b'inv_snl_paco', 59860 )
	ssc.data_set_number( data, b'inv_snl_pdco', 61130.78125 )
	ssc.data_set_number( data, b'inv_snl_pnt', 17.957999999999998 )
	ssc.data_set_number( data, b'inv_snl_pso', 97.213982000000001 )
	ssc.data_set_number( data, b'inv_snl_vdco', 630 )
	ssc.data_set_number( data, b'inv_snl_vdcmax', 800 )
	ssc.data_set_number( data, b'inv_cec_cg_c0', -3.1752000000000001e-06 )
	ssc.data_set_number( data, b'inv_cec_cg_c1', -5.1231399999999999e-05 )
	ssc.data_set_number( data, b'inv_cec_cg_c2', 0.00098359599999999999 )
	ssc.data_set_number( data, b'inv_cec_cg_c3', -0.0015077999999999999 )
	ssc.data_set_number( data, b'inv_cec_cg_paco', 3800 )
	ssc.data_set_number( data, b'inv_cec_cg_pdco', 3928.1100000000001 )
	ssc.data_set_number( data, b'inv_cec_cg_pnt', 0.98999999999999999 )
	ssc.data_set_number( data, b'inv_cec_cg_psco', 19.448399999999999 )
	ssc.data_set_number( data, b'inv_cec_cg_vdco', 398.49700000000001 )
	ssc.data_set_number( data, b'inv_cec_cg_vdcmax', 600 )
	ssc.data_set_number( data, b'inv_ds_paco', 4000 )
	ssc.data_set_number( data, b'inv_ds_eff', 96 )
	ssc.data_set_number( data, b'inv_ds_pnt', 1 )
	ssc.data_set_number( data, b'inv_ds_pso', 0 )
	ssc.data_set_number( data, b'inv_ds_vdco', 310 )
	ssc.data_set_number( data, b'inv_ds_vdcmax', 600 )
	ssc.data_set_number( data, b'inv_pd_paco', 4000 )
	ssc.data_set_number( data, b'inv_pd_pdco', 4210.5263671875 )
	inv_pd_partload =[ 0, 0.40400000000000003, 0.80800000000000005, 1.212, 1.6160000000000001, 2.02, 2.4239999999999999, 2.8279999999999998, 3.2320000000000002, 3.6360000000000001, 4.04, 4.444, 4.8479999999999999, 5.2519999999999998, 5.6559999999999997, 6.0599999999999996, 6.4640000000000004, 6.8680000000000003, 7.2720000000000002, 7.6760000000000002, 8.0800000000000001, 8.484, 8.8879999999999999, 9.2919999999999998, 9.6959999999999997, 10.1, 10.504, 10.907999999999999, 11.311999999999999, 11.715999999999999, 12.119999999999999, 12.523999999999999, 12.928000000000001, 13.332000000000001, 13.736000000000001, 14.140000000000001, 14.544, 14.948, 15.352, 15.756, 16.16, 16.564, 16.968, 17.372, 17.776, 18.18, 18.584, 18.988, 19.391999999999999, 19.795999999999999, 20.199999999999999, 20.603999999999999, 21.007999999999999, 21.411999999999999, 21.815999999999999, 22.219999999999999, 22.623999999999999, 23.027999999999999, 23.431999999999999, 23.835999999999999, 24.239999999999998, 24.643999999999998, 25.047999999999998, 25.452000000000002, 25.856000000000002, 26.260000000000002, 26.664000000000001, 27.068000000000001, 27.472000000000001, 27.876000000000001, 28.280000000000001, 28.684000000000001, 29.088000000000001, 29.492000000000001, 29.896000000000001, 30.300000000000001, 30.704000000000001, 31.108000000000001, 31.512, 31.916, 32.32, 32.723999999999997, 33.128, 33.531999999999996, 33.936, 34.340000000000003, 34.744, 35.148000000000003, 35.552, 35.956000000000003, 36.359999999999999, 36.764000000000003, 37.167999999999999, 37.572000000000003, 37.975999999999999, 38.380000000000003, 38.783999999999999, 39.188000000000002, 39.591999999999999, 39.996000000000002, 40.399999999999999, 40.804000000000002, 41.207999999999998, 41.612000000000002, 42.015999999999998, 42.420000000000002, 42.823999999999998, 43.228000000000002, 43.631999999999998, 44.036000000000001, 44.439999999999998, 44.844000000000001, 45.247999999999998, 45.652000000000001, 46.055999999999997, 46.460000000000001, 46.863999999999997, 47.268000000000001, 47.671999999999997, 48.076000000000001, 48.479999999999997, 48.884, 49.287999999999997, 49.692, 50.095999999999997, 50.5, 50.904000000000003, 51.308, 51.712000000000003, 52.116, 52.520000000000003, 52.923999999999999, 53.328000000000003, 53.731999999999999, 54.136000000000003, 54.539999999999999, 54.944000000000003, 55.347999999999999, 55.752000000000002, 56.155999999999999, 56.560000000000002, 56.963999999999999, 57.368000000000002, 57.771999999999998, 58.176000000000002, 58.579999999999998, 58.984000000000002, 59.387999999999998, 59.792000000000002, 60.195999999999998, 60.600000000000001, 61.003999999999998, 61.408000000000001, 61.811999999999998, 62.216000000000001, 62.619999999999997, 63.024000000000001, 63.427999999999997, 63.832000000000001, 64.236000000000004, 64.640000000000001, 65.043999999999997, 65.447999999999993, 65.852000000000004, 66.256, 66.659999999999997, 67.063999999999993, 67.468000000000004, 67.872, 68.275999999999996, 68.680000000000007, 69.084000000000003, 69.488, 69.891999999999996, 70.296000000000006, 70.700000000000003, 71.103999999999999, 71.507999999999996, 71.912000000000006, 72.316000000000003, 72.719999999999999, 73.123999999999995, 73.528000000000006, 73.932000000000002, 74.335999999999999, 74.739999999999995, 75.144000000000005, 75.548000000000002, 75.951999999999998, 76.355999999999995, 76.760000000000005, 77.164000000000001, 77.567999999999998, 77.971999999999994, 78.376000000000005, 78.780000000000001, 79.183999999999997, 79.587999999999994, 79.992000000000004, 80.396000000000001, 80.799999999999997, 81.203999999999994, 81.608000000000004, 82.012, 82.415999999999997, 82.819999999999993, 83.224000000000004, 83.628, 84.031999999999996, 84.436000000000007, 84.840000000000003, 85.244, 85.647999999999996, 86.052000000000007, 86.456000000000003, 86.859999999999999, 87.263999999999996, 87.668000000000006, 88.072000000000003, 88.475999999999999, 88.879999999999995, 89.284000000000006, 89.688000000000002, 90.091999999999999, 90.495999999999995, 90.900000000000006, 91.304000000000002, 91.707999999999998, 92.111999999999995, 92.516000000000005, 92.920000000000002, 93.323999999999998, 93.727999999999994, 94.132000000000005, 94.536000000000001, 94.939999999999998, 95.343999999999994, 95.748000000000005, 96.152000000000001, 96.555999999999997, 96.959999999999994, 97.364000000000004, 97.768000000000001, 98.171999999999997, 98.575999999999993, 98.980000000000004, 99.384, 99.787999999999997, 100.19199999999999, 100.596, 101 ];
	ssc.data_set_array( data, b'inv_pd_partload',  inv_pd_partload);
	inv_pd_efficiency =[ 0, 0, 34.420000000000002, 55.200000000000003, 65.590000000000003, 71.819999999999993, 75.969999999999999, 78.939999999999998, 81.170000000000002, 82.900000000000006, 84.280000000000001, 85.420000000000002, 86.359999999999999, 87.159999999999997, 87.840000000000003, 88.439999999999998, 88.950000000000003, 89.409999999999997, 89.819999999999993, 90.180000000000007, 90.510000000000005, 90.810000000000002, 91.079999999999998, 91.319999999999993, 91.549999999999997, 91.75, 91.950000000000003, 92.120000000000005, 92.290000000000006, 92.439999999999998, 92.579999999999998, 92.719999999999999, 92.840000000000003, 92.959999999999994, 93.069999999999993, 93.170000000000002, 93.269999999999996, 93.370000000000005, 93.450000000000003, 93.540000000000006, 93.620000000000005, 93.689999999999998, 93.760000000000005, 93.829999999999998, 93.900000000000006, 93.959999999999994, 94.019999999999996, 94.079999999999998, 94.129999999999995, 94.180000000000007, 94.230000000000004, 94.280000000000001, 94.329999999999998, 94.370000000000005, 94.420000000000002, 94.459999999999994, 94.5, 94.540000000000006, 94.569999999999993, 94.609999999999999, 94.640000000000001, 94.680000000000007, 94.709999999999994, 94.739999999999995, 94.769999999999996, 94.799999999999997, 94.829999999999998, 94.859999999999999, 94.890000000000001, 94.909999999999997, 94.939999999999998, 94.959999999999994, 94.980000000000004, 95.010000000000005, 95.030000000000001, 95.049999999999997, 95.069999999999993, 95.090000000000003, 95.109999999999999, 95.129999999999995, 95.150000000000006, 95.170000000000002, 95.189999999999998, 95.209999999999994, 95.230000000000004, 95.239999999999995, 95.260000000000005, 95.280000000000001, 95.290000000000006, 95.310000000000002, 95.319999999999993, 95.340000000000003, 95.349999999999994, 95.359999999999999, 95.379999999999995, 95.390000000000001, 95.400000000000006, 95.420000000000002, 95.430000000000007, 95.439999999999998, 95.450000000000003, 95.469999999999999, 95.480000000000004, 95.489999999999995, 95.5, 95.510000000000005, 95.519999999999996, 95.530000000000001, 95.540000000000006, 95.549999999999997, 95.560000000000002, 95.569999999999993, 95.579999999999998, 95.590000000000003, 95.599999999999994, 95.609999999999999, 95.620000000000005, 95.629999999999995, 95.640000000000001, 95.640000000000001, 95.650000000000006, 95.659999999999997, 95.670000000000002, 95.680000000000007, 95.680000000000007, 95.689999999999998, 95.700000000000003, 95.709999999999994, 95.709999999999994, 95.719999999999999, 95.730000000000004, 95.730000000000004, 95.739999999999995, 95.75, 95.75, 95.760000000000005, 95.769999999999996, 95.769999999999996, 95.780000000000001, 95.780000000000001, 95.790000000000006, 95.799999999999997, 95.799999999999997, 95.810000000000002, 95.810000000000002, 95.819999999999993, 95.819999999999993, 95.829999999999998, 95.829999999999998, 95.840000000000003, 95.840000000000003, 95.849999999999994, 95.849999999999994, 95.859999999999999, 95.859999999999999, 95.870000000000005, 95.870000000000005, 95.879999999999995, 95.879999999999995, 95.890000000000001, 95.890000000000001, 95.890000000000001, 95.900000000000006, 95.900000000000006, 95.909999999999997, 95.909999999999997, 95.909999999999997, 95.920000000000002, 95.920000000000002, 95.930000000000007, 95.930000000000007, 95.930000000000007, 95.939999999999998, 95.939999999999998, 95.939999999999998, 95.950000000000003, 95.950000000000003, 95.959999999999994, 95.959999999999994, 95.959999999999994, 95.969999999999999, 95.969999999999999, 95.969999999999999, 95.980000000000004, 95.980000000000004, 95.980000000000004, 95.980000000000004, 95.989999999999995, 95.989999999999995, 95.989999999999995, 96, 96, 96, 96.010000000000005, 96.010000000000005, 96.010000000000005, 96.010000000000005, 96.019999999999996, 96.019999999999996, 96.019999999999996, 96.019999999999996, 96.030000000000001, 96.030000000000001, 96.030000000000001, 96.030000000000001, 96.040000000000006, 96.040000000000006, 96.040000000000006, 96.040000000000006, 96.049999999999997, 96.049999999999997, 96.049999999999997, 96.049999999999997, 96.060000000000002, 96.060000000000002, 96.060000000000002, 96.060000000000002, 96.060000000000002, 96.069999999999993, 96.069999999999993, 96.069999999999993, 96.069999999999993, 96.069999999999993, 96.079999999999998, 96.079999999999998, 96.079999999999998, 96.079999999999998, 96.079999999999998, 96.090000000000003, 96.090000000000003, 96.090000000000003, 96.090000000000003, 96.090000000000003, 96.090000000000003, 96.099999999999994, 96.099999999999994, 96.099999999999994, 96.099999999999994, 96.099999999999994, 96.099999999999994, 96.109999999999999, 96.109999999999999, 96.109999999999999, 96.109999999999999, 96.109999999999999, 96.109999999999999, 96.120000000000005, 96.120000000000005, 96.120000000000005, 96.120000000000005, 96.120000000000005 ];
	ssc.data_set_array( data, b'inv_pd_efficiency',  inv_pd_efficiency);
	ssc.data_set_number( data, b'inv_pd_pnt', 0 )
	ssc.data_set_number( data, b'inv_pd_vdco', 310 )
	ssc.data_set_number( data, b'inv_pd_vdcmax', 600 )
	inv_tdc_cec_db = [[ 800,   28,   -0.02,   56,   0 ], [ 600,   52,   -0.037499999999999999,   60,   0 ], [ 390,   38,   -0.012500000000000001,   50,   -0.025000000000000001 ]];
	ssc.data_set_matrix( data, b'inv_tdc_cec_db', inv_tdc_cec_db );
	inv_tdc_cec_cg = [[ 800,   28,   -0.02,   56,   0 ], [ 600,   52,   -0.037499999999999999,   60,   0 ], [ 390,   38,   -0.012500000000000001,   50,   -0.025000000000000001 ]];
	ssc.data_set_matrix( data, b'inv_tdc_cec_cg', inv_tdc_cec_cg );
	inv_tdc_ds = [[ 800,   28,   -0.02,   56,   0 ], [ 600,   52,   -0.037499999999999999,   60,   0 ], [ 390,   38,   -0.012500000000000001,   50,   -0.025000000000000001 ]];
	ssc.data_set_matrix( data, b'inv_tdc_ds', inv_tdc_ds );
	inv_tdc_plc = [[ 800,   28,   -0.02,   56,   0 ], [ 600,   52,   -0.037499999999999999,   60,   0 ], [ 390,   38,   -0.012500000000000001,   50,   -0.025000000000000001 ]];
	ssc.data_set_matrix( data, b'inv_tdc_plc', inv_tdc_plc );
	ssc.data_set_number( data, b'en_batt', 1 )
	ssc.data_set_array_from_csv( data, b'load', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/load.csv');
	ssc.data_set_array_from_csv( data, b'crit_load', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/crit_load.csv');
	load_escalation =[ 0 ];
	ssc.data_set_array( data, b'load_escalation',  load_escalation);
	ssc.data_set_number( data, b'adjust:constant', 0 )
	ssc.data_set_number( data, b'dc_adjust:constant', 0 )
	ssc.data_set_number( data, b'batt_chem', 1 )
	ssc.data_set_number( data, b'inv_snl_eff_cec', 98.228355407714844 )
	ssc.data_set_number( data, b'inv_pd_eff', 95 )
	ssc.data_set_number( data, b'inv_cec_cg_eff_cec', 96.636390686035156 )
	ssc.data_set_number( data, b'batt_ac_or_dc', 1 )
	ssc.data_set_number( data, b'batt_dc_dc_efficiency', 99 )
	ssc.data_set_number( data, b'batt_dc_ac_efficiency', 96 )
	ssc.data_set_number( data, b'batt_ac_dc_efficiency', 96 )
	ssc.data_set_number( data, b'batt_meter_position', 0 )
	ssc.data_set_number( data, b'batt_inverter_efficiency_cutoff', 90 )
	batt_losses =[ 0 ];
	ssc.data_set_array( data, b'batt_losses',  batt_losses);
	batt_losses_charging =[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'batt_losses_charging',  batt_losses_charging);
	batt_losses_discharging =[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'batt_losses_discharging',  batt_losses_discharging);
	batt_losses_idle =[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'batt_losses_idle',  batt_losses_idle);
	ssc.data_set_number( data, b'batt_loss_choice', 0 )
	ssc.data_set_number( data, b'batt_current_choice', 1 )
	ssc.data_set_number( data, b'batt_computed_strings', (battcapacity/50 * 31) )
	ssc.data_set_number( data, b'batt_computed_series', 139 )
	ssc.data_set_number( data, b'batt_computed_bank_capacity', battcapacity )
	ssc.data_set_number( data, b'batt_current_charge_max', 200 )
	ssc.data_set_number( data, b'batt_current_discharge_max', 200 )
	ssc.data_set_number( data, b'batt_power_charge_max_kwdc', 100.08000183105469 )
	ssc.data_set_number( data, b'batt_power_discharge_max_kwdc', 100.08000183105469 )
	ssc.data_set_number( data, b'batt_power_charge_max_kwac', 104.25 )
	ssc.data_set_number( data, b'batt_power_discharge_max_kwac', 96.076797485351562 )
	ssc.data_set_number( data, b'batt_voltage_choice', 0 )
	ssc.data_set_number( data, b'batt_Vfull', 4.2000000000000002 )
	ssc.data_set_number( data, b'batt_Vexp', 3.5299999999999998 )
	ssc.data_set_number( data, b'batt_Vnom', 3.3420000000000001 )
	ssc.data_set_number( data, b'batt_Vnom_default', 3.6000000000000001 )
	ssc.data_set_number( data, b'batt_Qfull', 3.2000000000000002 )
	ssc.data_set_number( data, b'batt_Qfull_flow', 1200 )
	ssc.data_set_number( data, b'batt_Qexp', 2.5840001106262207 )
	ssc.data_set_number( data, b'batt_Qnom', 3.125999927520752 )
	ssc.data_set_number( data, b'batt_C_rate', 0.20000000000000001 )
	ssc.data_set_number( data, b'batt_resistance', 0.001155 )
	batt_voltage_matrix = [[ 0,   0 ]];
	ssc.data_set_matrix( data, b'batt_voltage_matrix', batt_voltage_matrix );
	ssc.data_set_number( data, b'LeadAcid_q20_computed', 1200 )
	ssc.data_set_number( data, b'LeadAcid_q10_computed', 1116 )
	ssc.data_set_number( data, b'LeadAcid_qn_computed', 720 )
	ssc.data_set_number( data, b'LeadAcid_tn', 1 )
	ssc.data_set_number( data, b'batt_initial_SOC', 50 )
	ssc.data_set_number( data, b'batt_minimum_SOC', 15 )
	ssc.data_set_number( data, b'batt_maximum_SOC', 95 )
	ssc.data_set_number( data, b'batt_minimum_modetime', 10 )
	batt_lifetime_matrix = [[ 10,   0,   100.85299999999999 ], [ 10,   1250,   94.884 ], [ 10,   2500,   88.914699999999996 ], [ 10,   3750,   82.945400000000006 ], [ 10,   5000,   76.976100000000002 ], [ 20,   0,   100.85299999999999 ], [ 20,   1250,   94.879800000000003 ], [ 20,   2500,   88.906300000000002 ], [ 20,   3750,   82.9328 ], [ 20,   5000,   76.959299999999999 ], [ 40,   0,   100.85299999999999 ], [ 40,   1250,   94.782200000000003 ], [ 40,   2500,   88.710999999999999 ], [ 40,   3750,   82.639700000000005 ], [ 40,   5000,   76.568200000000004 ], [ 80,   0,   100.85299999999999 ], [ 80,   1250,   92.483800000000002 ], [ 80,   2500,   84.054299999999998 ], [ 80,   3750,   75.560000000000002 ], [ 80,   5000,   66.995599999999996 ], [ 100,   0,   100.85299999999999 ], [ 100,   1250,   88.125600000000006 ], [ 100,   2500,   74.873199999999997 ], [ 100,   3750,   60.951099999999997 ], [ 100,   5000,   46.1312 ]];
	ssc.data_set_matrix( data, b'batt_lifetime_matrix', batt_lifetime_matrix );
	ssc.data_set_number( data, b'batt_calendar_choice', 1 )
	batt_calendar_lifetime_matrix = [[ 0,   100 ], [ 3650,   80 ], [ 7300,   50 ]];
	ssc.data_set_matrix( data, b'batt_calendar_lifetime_matrix', batt_calendar_lifetime_matrix );
	ssc.data_set_number( data, b'batt_calendar_q0', 1.02 )
	ssc.data_set_number( data, b'batt_calendar_a', 0.00266 )
	ssc.data_set_number( data, b'batt_calendar_b', -7280 )
	ssc.data_set_number( data, b'batt_calendar_c', 930 )
	ssc.data_set_number( data, b'batt_replacement_capacity', 50 )
	ssc.data_set_number( data, b'batt_replacement_option', 1 )
	batt_replacement_schedule_percent =[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'batt_replacement_schedule_percent',  batt_replacement_schedule_percent);
	om_replacement_cost1 =[ 206.66999999999999 ];
	ssc.data_set_array( data, b'om_replacement_cost1',  om_replacement_cost1);
	ssc.data_set_number( data, b'batt_mass', 5945.34619140625 )
	ssc.data_set_number( data, b'batt_surface_area', 25.14306640625 )
	ssc.data_set_number( data, b'batt_Cp', 1500 )
	ssc.data_set_number( data, b'batt_h_to_ambient', 7.5 )
	ssc.data_set_array_from_csv( data, b'batt_room_temperature_celsius', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/batt_room_temperature_celsius.csv');
	cap_vs_temp = [[ 0,   80.200000000000003 ], [ 23,   100 ], [ 30,   103.09999999999999 ], [ 45,   105.40000000000001 ]];
	ssc.data_set_matrix( data, b'cap_vs_temp', cap_vs_temp );
	dispatch_manual_charge =[ 1, 0, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'dispatch_manual_charge',  dispatch_manual_charge);
	dispatch_manual_discharge =[ 0, 1, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'dispatch_manual_discharge',  dispatch_manual_discharge);
	dispatch_manual_gridcharge =[ 0, 0, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'dispatch_manual_gridcharge',  dispatch_manual_gridcharge);
	dispatch_manual_percent_discharge =[ 25, 0 ];
	ssc.data_set_array( data, b'dispatch_manual_percent_discharge',  dispatch_manual_percent_discharge);
	dispatch_manual_percent_gridcharge =[ 0 ];
	ssc.data_set_array( data, b'dispatch_manual_percent_gridcharge',  dispatch_manual_percent_gridcharge);
	dispatch_manual_sched = [[ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3,   3,   3,   3,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3,   3,   3,   3,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3,   3,   3,   3,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3,   3,   3,   3,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   4,   4,   4,   4,   4,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   4,   4,   4,   4,   4,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   4,   4,   4,   4,   4,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   4,   4,   4,   4,   4,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   4,   4,   4,   4,   4,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   4,   4,   4,   4,   4,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3,   3,   3,   3,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3,   3,   3,   3,   1,   1,   1,   1 ]];
	ssc.data_set_matrix( data, b'dispatch_manual_sched', dispatch_manual_sched );
	dispatch_manual_sched_weekend = [[ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ]];
	ssc.data_set_matrix( data, b'dispatch_manual_sched_weekend', dispatch_manual_sched_weekend );
	ssc.data_set_array_from_csv( data, b'batt_target_power', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/batt_target_power.csv');
	batt_target_power_monthly =[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ];
	ssc.data_set_array( data, b'batt_target_power_monthly',  batt_target_power_monthly);
	ssc.data_set_number( data, b'batt_target_choice', 0 )
	batt_custom_dispatch =[ 0 ];
	ssc.data_set_array( data, b'batt_custom_dispatch',  batt_custom_dispatch);
	ssc.data_set_number( data, b'batt_dispatch_choice', 0 )
	ssc.data_set_number( data, b'batt_dispatch_auto_can_fuelcellcharge', 0 )
	ssc.data_set_number( data, b'batt_dispatch_auto_can_gridcharge', 0 )
	ssc.data_set_number( data, b'batt_dispatch_auto_can_charge', 1 )
	ssc.data_set_number( data, b'batt_cycle_cost_choice', 0 )
	batt_cycle_cost =[ 0.10000000000000001 ];
	ssc.data_set_array( data, b'batt_cycle_cost',  batt_cycle_cost);
	ssc.data_set_number( data, b'inflation_rate', 2.5 )
	rate_escalation =[ 0 ];
	ssc.data_set_array( data, b'rate_escalation',  rate_escalation);
	ssc.data_set_number( data, b'ur_metering_option', 0 )
	ssc.data_set_number( data, b'ur_nm_yearend_sell_rate', 0 )
	ssc.data_set_number( data, b'ur_nm_credit_month', 11 )
	ssc.data_set_number( data, b'ur_nm_credit_rollover', 0 )
	ssc.data_set_number( data, b'ur_monthly_fixed_charge', 30 )
	ssc.data_set_number( data, b'ur_monthly_min_charge', 0 )
	ssc.data_set_number( data, b'ur_annual_min_charge', 0 )
	ssc.data_set_number( data, b'ur_en_ts_sell_rate', 0 )
	ur_ts_sell_rate =[ 0 ];
	ssc.data_set_array( data, b'ur_ts_sell_rate',  ur_ts_sell_rate);
	ssc.data_set_number( data, b'ur_en_ts_buy_rate', 0 )
	ssc.data_set_array_from_csv( data, b'ur_ts_buy_rate', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/ur_ts_buy_rate.csv');
	ur_ec_sched_weekday = [[ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   3,   3,   3,   3,   3,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   3,   3,   3,   3,   3,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   3,   3,   3,   3,   3,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   3,   3,   3,   3,   3,   4,   4,   4,   4 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   3,   3,   3,   3,   3,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   3,   3,   3,   3,   3,   4,   4,   4,   4 ]];
	ssc.data_set_matrix( data, b'ur_ec_sched_weekday', ur_ec_sched_weekday );
	ur_ec_sched_weekend = [[ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4 ], [ 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4 ]];
	ssc.data_set_matrix( data, b'ur_ec_sched_weekend', ur_ec_sched_weekend );
	ur_ec_tou_mat = [[ 1,   1,   9.9999999999999998e+37,   0,   0.050000000000000003,   0 ], [ 2,   1,   9.9999999999999998e+37,   0,   0.074999999999999997,   0 ], [ 3,   1,   9.9999999999999998e+37,   0,   0.059999999999999998,   0 ], [ 4,   1,   9.9999999999999998e+37,   0,   0.050000000000000003,   0 ]];
	ssc.data_set_matrix( data, b'ur_ec_tou_mat', ur_ec_tou_mat );
	ssc.data_set_number( data, b'ur_dc_enable', 1 )
	ur_dc_sched_weekday = [[ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   2,   2,   2,   2 ]];
	ssc.data_set_matrix( data, b'ur_dc_sched_weekday', ur_dc_sched_weekday );
	ur_dc_sched_weekend = [[ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ], [ 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2 ]];
	ssc.data_set_matrix( data, b'ur_dc_sched_weekend', ur_dc_sched_weekend );
	ur_dc_tou_mat = [[ 1,   1,   100,   25 ], [ 1,   2,   9.9999999999999998e+37,   35 ], [ 2,   1,   100,   25 ], [ 2,   2,   9.9999999999999998e+37,   15 ]];
	ssc.data_set_matrix( data, b'ur_dc_tou_mat', ur_dc_tou_mat );
	ur_dc_flat_mat = [[ 0,   1,   9.9999999999999998e+37,   0 ], [ 1,   1,   9.9999999999999998e+37,   0 ], [ 2,   1,   9.9999999999999998e+37,   0 ], [ 3,   1,   9.9999999999999998e+37,   0 ], [ 4,   1,   9.9999999999999998e+37,   0 ], [ 5,   1,   9.9999999999999998e+37,   0 ], [ 6,   1,   9.9999999999999998e+37,   0 ], [ 7,   1,   9.9999999999999998e+37,   0 ], [ 8,   1,   9.9999999999999998e+37,   0 ], [ 9,   1,   9.9999999999999998e+37,   0 ], [ 10,   1,   9.9999999999999998e+37,   0 ], [ 11,   1,   9.9999999999999998e+37,   0 ]];
	ssc.data_set_matrix( data, b'ur_dc_flat_mat', ur_dc_flat_mat );
	module = ssc.module_create(b'pvsamv1')	
	ssc.module_exec_set_print( 0 );
	if ssc.module_exec(module, data) == 0:
		print ('pvsamv1 simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" );
	ssc.module_free(module)
	ssc.data_set_number( data, b'enable_interconnection_limit', 0 )
	ssc.data_set_number( data, b'grid_interconnection_limit_kwac', 100000 )
	ssc.data_set_array_from_csv( data, b'grid_curtailment', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/BatterySweep/grid_curtailment.csv');
	module = ssc.module_create(b'grid')	
	ssc.module_exec_set_print( 0 );
	if ssc.module_exec(module, data) == 0:
		print ('grid simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" );
	ssc.module_free(module)
	degradation =[ 0 ];
	ssc.data_set_array( data, b'degradation',  degradation);
	module = ssc.module_create(b'utilityrate5')	
	ssc.module_exec_set_print( 0 );
	if ssc.module_exec(module, data) == 0:
		print ('utilityrate5 simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" );
	ssc.module_free(module)
	federal_tax_rate =[ 21 ];
	ssc.data_set_array( data, b'federal_tax_rate',  federal_tax_rate);
	state_tax_rate =[ 7 ];
	ssc.data_set_array( data, b'state_tax_rate',  state_tax_rate);
	ssc.data_set_number( data, b'property_tax_rate', 0 )
	ssc.data_set_number( data, b'prop_tax_cost_assessed_percent', 100 )
	ssc.data_set_number( data, b'prop_tax_assessed_decline', 0 )
	ssc.data_set_number( data, b'real_discount_rate', 6.4000000000000004 )
	ssc.data_set_number( data, b'insurance_rate', 0 )
	ssc.data_set_number( data, b'loan_term', 25 )
	ssc.data_set_number( data, b'loan_rate', 5 )
	ssc.data_set_number( data, b'debt_fraction', 60 )
	om_fixed =[ 0 ];
	ssc.data_set_array( data, b'om_fixed',  om_fixed);
	ssc.data_set_number( data, b'om_fixed_escal', 0 )
	om_production =[ 0 ];
	ssc.data_set_array( data, b'om_production',  om_production);
	ssc.data_set_number( data, b'om_production_escal', 0 )
	om_capacity =[ 19 ];
	ssc.data_set_array( data, b'om_capacity',  om_capacity);
	ssc.data_set_number( data, b'om_capacity_escal', 0 )
	om_fuel_cost =[ 0 ];
	ssc.data_set_array( data, b'om_fuel_cost',  om_fuel_cost);
	ssc.data_set_number( data, b'om_fuel_cost_escal', 0 )
	ssc.data_set_number( data, b'om_replacement_cost_escal', 0 )
	ssc.data_set_number( data, b'depr_fed_type', 1 )
	ssc.data_set_number( data, b'depr_fed_sl_years', 7 )
	depr_fed_custom =[ 0 ];
	ssc.data_set_array( data, b'depr_fed_custom',  depr_fed_custom);
	ssc.data_set_number( data, b'depr_sta_type', 1 )
	ssc.data_set_number( data, b'depr_sta_sl_years', 7 )
	depr_sta_custom =[ 0 ];
	ssc.data_set_array( data, b'depr_sta_custom',  depr_sta_custom);
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
	ptc_fed_amount =[ 0 ];
	ssc.data_set_array( data, b'ptc_fed_amount',  ptc_fed_amount);
	ssc.data_set_number( data, b'ptc_fed_term', 10 )
	ssc.data_set_number( data, b'ptc_fed_escal', 0 )
	ptc_sta_amount =[ 0 ];
	ssc.data_set_array( data, b'ptc_sta_amount',  ptc_sta_amount);
	ssc.data_set_number( data, b'ptc_sta_term', 10 )
	ssc.data_set_number( data, b'ptc_sta_escal', 0 )
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
	pbi_fed_amount =[ 0 ];
	ssc.data_set_array( data, b'pbi_fed_amount',  pbi_fed_amount);
	ssc.data_set_number( data, b'pbi_fed_term', 0 )
	ssc.data_set_number( data, b'pbi_fed_escal', 0 )
	ssc.data_set_number( data, b'pbi_fed_tax_fed', 1 )
	ssc.data_set_number( data, b'pbi_fed_tax_sta', 1 )
	pbi_sta_amount =[ 0 ];
	ssc.data_set_array( data, b'pbi_sta_amount',  pbi_sta_amount);
	ssc.data_set_number( data, b'pbi_sta_term', 0 )
	ssc.data_set_number( data, b'pbi_sta_escal', 0 )
	ssc.data_set_number( data, b'pbi_sta_tax_fed', 1 )
	ssc.data_set_number( data, b'pbi_sta_tax_sta', 1 )
	pbi_uti_amount =[ 0 ];
	ssc.data_set_array( data, b'pbi_uti_amount',  pbi_uti_amount);
	ssc.data_set_number( data, b'pbi_uti_term', 0 )
	ssc.data_set_number( data, b'pbi_uti_escal', 0 )
	ssc.data_set_number( data, b'pbi_uti_tax_fed', 1 )
	ssc.data_set_number( data, b'pbi_uti_tax_sta', 1 )
	pbi_oth_amount =[ 0 ];
	ssc.data_set_array( data, b'pbi_oth_amount',  pbi_oth_amount);
	ssc.data_set_number( data, b'pbi_oth_term', 0 )
	ssc.data_set_number( data, b'pbi_oth_escal', 0 )
	ssc.data_set_number( data, b'pbi_oth_tax_fed', 1 )
	ssc.data_set_number( data, b'pbi_oth_tax_sta', 1 )
	ssc.data_set_number( data, b'battery_per_kWh', 206.66999999999999 )
	ssc.data_set_number( data, b'total_installed_cost', 1029040.0625 )
	ssc.data_set_number( data, b'salvage_percentage', 0 )
	module = ssc.module_create(b'cashloan')	
	ssc.module_exec_set_print( 0 );
	if ssc.module_exec(module, data) == 0:
		print ('cashloan simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" );
	ssc.module_free(module)
# 	annual_energy = ssc.data_get_number(data, b'annual_energy');
# 	print ('Annual energy (year 1) = ', annual_energy)
# 	capacity_factor = ssc.data_get_number(data, b'capacity_factor');
# 	print ('Capacity factor (year 1) = ', capacity_factor)
# 	kwh_per_kw = ssc.data_get_number(data, b'kwh_per_kw');
# 	print ('Energy yield (year 1) = ', kwh_per_kw)
# 	performance_ratio = ssc.data_get_number(data, b'performance_ratio');
# 	print ('Performance ratio (year 1) = ', performance_ratio)
# 	average_battery_roundtrip_efficiency = ssc.data_get_number(data, b'average_battery_roundtrip_efficiency');
# 	print ('Battery roundtrip efficiency = ', average_battery_roundtrip_efficiency)
# 	batt_system_charge_percent = ssc.data_get_number(data, b'batt_system_charge_percent');
# 	print ('Battery charge energy from system = ', batt_system_charge_percent)
# 	lcoe_nom = ssc.data_get_number(data, b'lcoe_nom');
# 	print ('Levelized COE (nominal) = ', lcoe_nom)
# 	lcoe_real = ssc.data_get_number(data, b'lcoe_real');
# 	print ('Levelized COE (real) = ', lcoe_real)
# 	elec_cost_without_system_year1 = ssc.data_get_number(data, b'elec_cost_without_system_year1');
# 	print ('Electricity bill without system (year 1) = ', elec_cost_without_system_year1)
# 	elec_cost_with_system_year1 = ssc.data_get_number(data, b'elec_cost_with_system_year1');
# 	print ('Electricity bill with system (year 1) = ', elec_cost_with_system_year1)
# 	savings_year1 = ssc.data_get_number(data, b'savings_year1');
# 	print ('Net savings with system (year 1) = ', savings_year1)
# 	npv = ssc.data_get_number(data, b'npv');
# 	print ('Net present value = ', npv)
	payback = ssc.data_get_number(data, b'payback');
	print ('Simple payback period = ', payback)
# 	discounted_payback = ssc.data_get_number(data, b'discounted_payback');
# 	print ('Discounted payback period = ', discounted_payback)
# 	adjusted_installed_cost = ssc.data_get_number(data, b'adjusted_installed_cost');
# 	print ('Net capital cost = ', adjusted_installed_cost)
# 	first_cost = ssc.data_get_number(data, b'first_cost');
# 	print ('Equity = ', first_cost)
# 	loan_amount = ssc.data_get_number(data, b'loan_amount');
# 	print ('Debt = ', loan_amount)
	ssc.data_free(data);
	return(payback)
    