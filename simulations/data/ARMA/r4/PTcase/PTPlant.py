from PySSC import PySSC

def ptrun(weather):
	ssc = PySSC()
	print ('Current folder = C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/PTcase')
	print ('SSC Version = ', ssc.version())
	print ('SSC Build Information = ', ssc.build_info().decode("utf - 8"))
	ssc.module_exec_set_print(0)
	data = ssc.data_create()
	ssc.data_set_string( data, b'file_name', weather );
	ssc.data_set_number( data, b'nSCA', 8 )
	ssc.data_set_number( data, b'nHCEt', 4 )
	ssc.data_set_number( data, b'nColt', 4 )
	ssc.data_set_number( data, b'nHCEVar', 4 )
	ssc.data_set_number( data, b'nLoops', 181 )
	ssc.data_set_number( data, b'FieldConfig', 2 )
	ssc.data_set_number( data, b'include_fixed_power_block_runner', 1 )
	ssc.data_set_number( data, b'L_power_block_piping', 50 )
	ssc.data_set_number( data, b'eta_pump', 0.84999999999999998 )
	ssc.data_set_number( data, b'Fluid', 21 )
	ssc.data_set_number( data, b'accept_loc', 1 )
	ssc.data_set_number( data, b'HDR_rough', 4.57e-05 )
	ssc.data_set_number( data, b'theta_stow', 170 )
	ssc.data_set_number( data, b'theta_dep', 10 )
	ssc.data_set_number( data, b'Row_Distance', 15 )
	ssc.data_set_number( data, b'T_loop_in_des', 293 )
	ssc.data_set_number( data, b'T_loop_out', 391 )
	ssc.data_set_number( data, b'm_dot_htfmin', 1 )
	ssc.data_set_number( data, b'm_dot_htfmax', 12 )
	field_fl_props = [[ 0 ]];
	ssc.data_set_matrix( data, b'field_fl_props', field_fl_props );
	ssc.data_set_number( data, b'T_fp', 150 )
	ssc.data_set_number( data, b'I_bn_des', 950 )
	ssc.data_set_number( data, b'Pipe_hl_coef', 0.45000000000000001 )
	ssc.data_set_number( data, b'SCA_drives_elec', 125 )
	ssc.data_set_number( data, b'tilt', 0 )
	ssc.data_set_number( data, b'azimuth', 0 )
	ssc.data_set_number( data, b'wind_stow_speed', 25 )
	ssc.data_set_number( data, b'accept_mode', 0 )
	ssc.data_set_number( data, b'accept_init', 0 )
	ssc.data_set_number( data, b'solar_mult', 2 )
	ssc.data_set_number( data, b'mc_bal_hot', 0.20000000000000001 )
	ssc.data_set_number( data, b'mc_bal_cold', 0.20000000000000001 )
	ssc.data_set_number( data, b'mc_bal_sca', 4.5 )
	W_aperture =[ 6, 6, 6, 6 ];
	ssc.data_set_array( data, b'W_aperture',  W_aperture);
	A_aperture =[ 656, 656, 656, 656 ];
	ssc.data_set_array( data, b'A_aperture',  A_aperture);
	TrackingError =[ 0.98799999999999999, 0.98799999999999999, 0.98799999999999999, 0.98799999999999999 ];
	ssc.data_set_array( data, b'TrackingError',  TrackingError);
	GeomEffects =[ 0.95199999999999996, 0.95199999999999996, 0.95199999999999996, 0.95199999999999996 ];
	ssc.data_set_array( data, b'GeomEffects',  GeomEffects);
	Rho_mirror_clean =[ 0.93000000000000005, 0.93000000000000005, 0.93000000000000005, 0.93000000000000005 ];
	ssc.data_set_array( data, b'Rho_mirror_clean',  Rho_mirror_clean);
	Dirt_mirror =[ 0.96999999999999997, 0.96999999999999997, 0.96999999999999997, 0.96999999999999997 ];
	ssc.data_set_array( data, b'Dirt_mirror',  Dirt_mirror);
	Error =[ 1, 1, 1, 1 ];
	ssc.data_set_array( data, b'Error',  Error);
	Ave_Focal_Length =[ 2.1499999999999999, 2.1499999999999999, 2.1499999999999999, 2.1499999999999999 ];
	ssc.data_set_array( data, b'Ave_Focal_Length',  Ave_Focal_Length);
	L_SCA =[ 115, 115, 115, 115 ];
	ssc.data_set_array( data, b'L_SCA',  L_SCA);
	L_aperture =[ 14.375, 14.375, 14.375, 14.375 ];
	ssc.data_set_array( data, b'L_aperture',  L_aperture);
	ColperSCA =[ 8, 8, 8, 8 ];
	ssc.data_set_array( data, b'ColperSCA',  ColperSCA);
	Distance_SCA =[ 1, 1, 1, 1 ];
	ssc.data_set_array( data, b'Distance_SCA',  Distance_SCA);
	IAM_matrix = [[ 1,   0.0327,   -0.1351 ], [ 1,   0.0327,   -0.1351 ], [ 1,   0.0327,   -0.1351 ], [ 1,   0.0327,   -0.1351 ]];
	ssc.data_set_matrix( data, b'IAM_matrix', IAM_matrix );
	HCE_FieldFrac = [[ 0.98499999999999999,   0.01,   0.0050000000000000001,   0 ], [ 1,   0,   0,   0 ], [ 1,   0,   0,   0 ], [ 1,   0,   0,   0 ]];
	ssc.data_set_matrix( data, b'HCE_FieldFrac', HCE_FieldFrac );
	D_2 = [[ 0.075999999999999998,   0.075999999999999998,   0.075999999999999998,   0.075999999999999998 ], [ 0.075999999999999998,   0.075999999999999998,   0.075999999999999998,   0.075999999999999998 ], [ 0.075999999999999998,   0.075999999999999998,   0.075999999999999998,   0.075999999999999998 ], [ 0.075999999999999998,   0.075999999999999998,   0.075999999999999998,   0.075999999999999998 ]];
	ssc.data_set_matrix( data, b'D_2', D_2 );
	D_3 = [[ 0.080000000000000002,   0.080000000000000002,   0.080000000000000002,   0.080000000000000002 ], [ 0.080000000000000002,   0.080000000000000002,   0.080000000000000002,   0.080000000000000002 ], [ 0.080000000000000002,   0.080000000000000002,   0.080000000000000002,   0.080000000000000002 ], [ 0.080000000000000002,   0.080000000000000002,   0.080000000000000002,   0.080000000000000002 ]];
	ssc.data_set_matrix( data, b'D_3', D_3 );
	D_4 = [[ 0.115,   0.115,   0.115,   0.115 ], [ 0.115,   0.115,   0.115,   0.115 ], [ 0.115,   0.115,   0.115,   0.115 ], [ 0.115,   0.115,   0.115,   0.115 ]];
	ssc.data_set_matrix( data, b'D_4', D_4 );
	D_5 = [[ 0.12,   0.12,   0.12,   0.12 ], [ 0.12,   0.12,   0.12,   0.12 ], [ 0.12,   0.12,   0.12,   0.12 ], [ 0.12,   0.12,   0.12,   0.12 ]];
	ssc.data_set_matrix( data, b'D_5', D_5 );
	D_p = [[ 0,   0,   0,   0 ], [ 0,   0,   0,   0 ], [ 0,   0,   0,   0 ], [ 0,   0,   0,   0 ]];
	ssc.data_set_matrix( data, b'D_p', D_p );
	Flow_type = [[ 1,   1,   1,   1 ], [ 1,   1,   1,   1 ], [ 1,   1,   1,   1 ], [ 1,   1,   1,   1 ]];
	ssc.data_set_matrix( data, b'Flow_type', Flow_type );
	Rough = [[ 4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05 ], [ 4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05 ], [ 4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05 ], [ 4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05,   4.5000000000000003e-05 ]];
	ssc.data_set_matrix( data, b'Rough', Rough );
	alpha_env = [[ 0.02,   0.02,   0,   0 ], [ 0.02,   0.02,   0,   0 ], [ 0.02,   0.02,   0,   0 ], [ 0.02,   0.02,   0,   0 ]];
	ssc.data_set_matrix( data, b'alpha_env', alpha_env );
	epsilon_3_11 = [[ 100,   0.064000000000000001 ], [ 150,   0.066500000000000004 ], [ 200,   0.070000000000000007 ], [ 250,   0.074499999999999997 ], [ 300,   0.080000000000000002 ], [ 350,   0.086499999999999994 ], [ 400,   0.094 ], [ 450,   0.10249999999999999 ], [ 500,   0.112 ]];
	ssc.data_set_matrix( data, b'epsilon_3_11', epsilon_3_11 );
	epsilon_3_12 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_12', epsilon_3_12 );
	epsilon_3_13 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_13', epsilon_3_13 );
	epsilon_3_14 = [[ 0 ]];
	ssc.data_set_matrix( data, b'epsilon_3_14', epsilon_3_14 );
	epsilon_3_21 = [[ 100,   0.064000000000000001 ], [ 150,   0.066500000000000004 ], [ 200,   0.070000000000000007 ], [ 250,   0.074499999999999997 ], [ 300,   0.080000000000000002 ], [ 350,   0.086499999999999994 ], [ 400,   0.094 ], [ 450,   0.10249999999999999 ], [ 500,   0.112 ]];
	ssc.data_set_matrix( data, b'epsilon_3_21', epsilon_3_21 );
	epsilon_3_22 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_22', epsilon_3_22 );
	epsilon_3_23 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_23', epsilon_3_23 );
	epsilon_3_24 = [[ 0 ]];
	ssc.data_set_matrix( data, b'epsilon_3_24', epsilon_3_24 );
	epsilon_3_31 = [[ 100,   0.064000000000000001 ], [ 150,   0.066500000000000004 ], [ 200,   0.070000000000000007 ], [ 250,   0.074499999999999997 ], [ 300,   0.080000000000000002 ], [ 350,   0.086499999999999994 ], [ 400,   0.094 ], [ 450,   0.10249999999999999 ], [ 500,   0.112 ]];
	ssc.data_set_matrix( data, b'epsilon_3_31', epsilon_3_31 );
	epsilon_3_32 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_32', epsilon_3_32 );
	epsilon_3_33 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_33', epsilon_3_33 );
	epsilon_3_34 = [[ 0 ]];
	ssc.data_set_matrix( data, b'epsilon_3_34', epsilon_3_34 );
	epsilon_3_41 = [[ 100,   0.064000000000000001 ], [ 150,   0.066500000000000004 ], [ 200,   0.070000000000000007 ], [ 250,   0.074499999999999997 ], [ 300,   0.080000000000000002 ], [ 350,   0.086499999999999994 ], [ 400,   0.094 ], [ 450,   0.10249999999999999 ], [ 500,   0.112 ]];
	ssc.data_set_matrix( data, b'epsilon_3_41', epsilon_3_41 );
	epsilon_3_42 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_42', epsilon_3_42 );
	epsilon_3_43 = [[ 0.65000000000000002 ]];
	ssc.data_set_matrix( data, b'epsilon_3_43', epsilon_3_43 );
	epsilon_3_44 = [[ 0 ]];
	ssc.data_set_matrix( data, b'epsilon_3_44', epsilon_3_44 );
	alpha_abs = [[ 0.96299999999999997,   0.96299999999999997,   0.80000000000000004,   0 ], [ 0.96299999999999997,   0.96299999999999997,   0.80000000000000004,   0 ], [ 0.96299999999999997,   0.96299999999999997,   0.80000000000000004,   0 ], [ 0.96299999999999997,   0.96299999999999997,   0.80000000000000004,   0 ]];
	ssc.data_set_matrix( data, b'alpha_abs', alpha_abs );
	Tau_envelope = [[ 0.96399999999999997,   0.96399999999999997,   1,   0 ], [ 0.96399999999999997,   0.96399999999999997,   1,   0 ], [ 0.96399999999999997,   0.96399999999999997,   1,   0 ], [ 0.96399999999999997,   0.96399999999999997,   1,   0 ]];
	ssc.data_set_matrix( data, b'Tau_envelope', Tau_envelope );
	EPSILON_4 = [[ 0.85999999999999999,   0.85999999999999999,   1,   0 ], [ 0.85999999999999999,   0.85999999999999999,   1,   0 ], [ 0.85999999999999999,   0.85999999999999999,   1,   0 ], [ 0.85999999999999999,   0.85999999999999999,   1,   0 ]];
	ssc.data_set_matrix( data, b'EPSILON_4', EPSILON_4 );
	EPSILON_5 = [[ 0.85999999999999999,   0.85999999999999999,   1,   0 ], [ 0.85999999999999999,   0.85999999999999999,   1,   0 ], [ 0.85999999999999999,   0.85999999999999999,   1,   0 ], [ 0.85999999999999999,   0.85999999999999999,   1,   0 ]];
	ssc.data_set_matrix( data, b'EPSILON_5', EPSILON_5 );
	GlazingIntactIn = [[ 1,   1,   0,   1 ], [ 1,   1,   0,   1 ], [ 1,   1,   0,   1 ], [ 1,   1,   0,   1 ]];
	ssc.data_set_matrix( data, b'GlazingIntactIn', GlazingIntactIn );
	P_a = [[ 0.0001,   750,   750,   0 ], [ 0.0001,   750,   750,   0 ], [ 0.0001,   750,   750,   0 ], [ 0.0001,   750,   750,   0 ]];
	ssc.data_set_matrix( data, b'P_a', P_a );
	AnnulusGas = [[ 27,   1,   1,   27 ], [ 27,   1,   1,   27 ], [ 27,   1,   1,   27 ], [ 27,   1,   1,   27 ]];
	ssc.data_set_matrix( data, b'AnnulusGas', AnnulusGas );
	AbsorberMaterial = [[ 1,   1,   1,   1 ], [ 1,   1,   1,   1 ], [ 1,   1,   1,   1 ], [ 1,   1,   1,   1 ]];
	ssc.data_set_matrix( data, b'AbsorberMaterial', AbsorberMaterial );
	Shadowing = [[ 0.93500000000000005,   0.93500000000000005,   0.93500000000000005,   0.96299999999999997 ], [ 0.93500000000000005,   0.93500000000000005,   0.93500000000000005,   0.96299999999999997 ], [ 0.93500000000000005,   0.93500000000000005,   0.93500000000000005,   0.96299999999999997 ], [ 0.93500000000000005,   0.93500000000000005,   0.93500000000000005,   0.96299999999999997 ]];
	ssc.data_set_matrix( data, b'Shadowing', Shadowing );
	Dirt_HCE = [[ 0.97999999999999998,   0.97999999999999998,   1,   0.97999999999999998 ], [ 0.97999999999999998,   0.97999999999999998,   1,   0.97999999999999998 ], [ 0.97999999999999998,   0.97999999999999998,   1,   0.97999999999999998 ], [ 0.97999999999999998,   0.97999999999999998,   1,   0.97999999999999998 ]];
	ssc.data_set_matrix( data, b'Dirt_HCE', Dirt_HCE );
	Design_loss = [[ 190,   1270,   1500,   0 ], [ 190,   1270,   1500,   0 ], [ 190,   1270,   1500,   0 ], [ 190,   1270,   1500,   0 ]];
	ssc.data_set_matrix( data, b'Design_loss', Design_loss );
	SCAInfoArray = [[ 1,   1 ], [ 1,   1 ], [ 1,   1 ], [ 1,   1 ], [ 1,   1 ], [ 1,   1 ], [ 1,   1 ], [ 1,   1 ]];
	ssc.data_set_matrix( data, b'SCAInfoArray', SCAInfoArray );
	SCADefocusArray =[ 8, 7, 6, 5, 4, 3, 2, 1 ];
	ssc.data_set_array( data, b'SCADefocusArray',  SCADefocusArray);
	ssc.data_set_number( data, b'rec_su_delay', 0.20000000000000001 )
	ssc.data_set_number( data, b'rec_qf_delay', 0.25 )
	ssc.data_set_number( data, b'p_start', 0.021000000000000001 )
	ssc.data_set_number( data, b'pc_config', 0 )
	ssc.data_set_number( data, b'P_ref', 111 )
	ssc.data_set_number( data, b'eta_ref', 0.35599999999999998 )
	ssc.data_set_number( data, b'cycle_max_frac', 1.05 )
	ssc.data_set_number( data, b'cycle_cutoff_frac', 0.20000000000000001 )
	ssc.data_set_number( data, b'q_sby_frac', 0.20000000000000001 )
	ssc.data_set_number( data, b'startup_time', 0.5 )
	ssc.data_set_number( data, b'startup_frac', 0.20000000000000001 )
	ssc.data_set_number( data, b'pb_pump_coef', 0.55000000000000004 )
	ssc.data_set_number( data, b'dT_cw_ref', 10 )
	ssc.data_set_number( data, b'T_amb_des', 42 )
	ssc.data_set_number( data, b'P_boil', 100 )
	ssc.data_set_number( data, b'CT', 2 )
	ssc.data_set_number( data, b'tech_type', 1 )
	ssc.data_set_number( data, b'T_approach', 5 )
	ssc.data_set_number( data, b'T_ITD_des', 16 )
	ssc.data_set_number( data, b'P_cond_ratio', 1.0027999999999999 )
	ssc.data_set_number( data, b'pb_bd_frac', 0.02 )
	ssc.data_set_number( data, b'P_cond_min', 1.25 )
	ssc.data_set_number( data, b'n_pl_inc', 8 )
	F_wc =[ 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
	ssc.data_set_array( data, b'F_wc',  F_wc);
	ssc.data_set_number( data, b'ud_f_W_dot_cool_des', 0 )
	ssc.data_set_number( data, b'ud_m_dot_water_cool_des', 0 )
	ssc.data_set_matrix_from_csv( data, b'ud_ind_od', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/PTcase/ud_ind_od.csv');
	ssc.data_set_number( data, b'store_fluid', 18 )
	store_fl_props = [[ 1 ]];
	ssc.data_set_matrix( data, b'store_fl_props', store_fl_props );
	ssc.data_set_number( data, b'is_hx', 1 )
	ssc.data_set_number( data, b'tshours', 6 )
	ssc.data_set_number( data, b'h_tank', 12 )
	ssc.data_set_number( data, b'u_tank', 0.40000000000000002 )
	ssc.data_set_number( data, b'tank_pairs', 1 )
	ssc.data_set_number( data, b'hot_tank_Thtr', 365 )
	ssc.data_set_number( data, b'hot_tank_max_heat', 25 )
	ssc.data_set_number( data, b'cold_tank_Thtr', 250 )
	ssc.data_set_number( data, b'cold_tank_max_heat', 25 )
	ssc.data_set_number( data, b'dt_hot', 5 )
	ssc.data_set_number( data, b'h_tank_min', 1 )
	ssc.data_set_number( data, b'init_hot_htf_percent', 30 )
	weekday_schedule = [[ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ]];
	ssc.data_set_matrix( data, b'weekday_schedule', weekday_schedule );
	weekend_schedule = [[ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ]];
	ssc.data_set_matrix( data, b'weekend_schedule', weekend_schedule );
	dispatch_sched_weekday = [[ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ]];
	ssc.data_set_matrix( data, b'dispatch_sched_weekday', dispatch_sched_weekday );
	dispatch_sched_weekend = [[ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 3,   3,   3,   3,   3,   3,   3,   3,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   2,   2,   2,   3,   3,   3 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ], [ 6,   6,   6,   6,   6,   6,   5,   5,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5 ]];
	ssc.data_set_matrix( data, b'dispatch_sched_weekend', dispatch_sched_weekend );
	ssc.data_set_number( data, b'is_tod_pc_target_also_pc_max', 0 )
	ssc.data_set_number( data, b'is_dispatch', 0 )
	ssc.data_set_number( data, b'disp_frequency', 24 )
	ssc.data_set_number( data, b'disp_horizon', 48 )
	ssc.data_set_number( data, b'disp_max_iter', 35000 )
	ssc.data_set_number( data, b'disp_timeout', 5 )
	ssc.data_set_number( data, b'disp_mip_gap', 0.001 )
	ssc.data_set_number( data, b'disp_time_weighting', 0.98999999999999999 )
	ssc.data_set_number( data, b'disp_rsu_cost', 950 )
	ssc.data_set_number( data, b'disp_csu_cost', 10000 )
	ssc.data_set_number( data, b'disp_pen_delta_w', 0.10000000000000001 )
	ssc.data_set_number( data, b'is_wlim_series', 0 )
	ssc.data_set_array_from_csv( data, b'wlim_series', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/PTcase/wlim_series.csv');
	f_turb_tou_periods =[ 1.05, 1, 1, 1, 1, 1, 1, 1, 1 ];
	ssc.data_set_array( data, b'f_turb_tou_periods',  f_turb_tou_periods);
	ssc.data_set_number( data, b'ppa_multiplier_model', 0 )
	ssc.data_set_array_from_csv( data, b'dispatch_factors_ts', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/PTcase/dispatch_factors_ts.csv');
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
	dispatch_series =[ 0 ];
	ssc.data_set_array( data, b'dispatch_series',  dispatch_series);
	ssc.data_set_number( data, b'pb_fixed_par', 0.0054999999999999997 )
	bop_array =[ 0, 1, 0.48299999999999998, 0.51700000000000002, 0 ];
	ssc.data_set_array( data, b'bop_array',  bop_array);
	aux_array =[ 0.02273, 1, 0.48299999999999998, 0.51700000000000002, 0 ];
	ssc.data_set_array( data, b'aux_array',  aux_array);
	ssc.data_set_number( data, b'gross_net_conversion_factor', 0.90000000000000002 )
	ssc.data_set_number( data, b'water_usage_per_wash', 0.69999999999999996 )
	ssc.data_set_number( data, b'washing_frequency', 63 )
	ssc.data_set_number( data, b'calc_design_pipe_vals', 1 )
	ssc.data_set_number( data, b'V_hdr_cold_max', 3 )
	ssc.data_set_number( data, b'V_hdr_cold_min', 2 )
	ssc.data_set_number( data, b'V_hdr_hot_max', 3 )
	ssc.data_set_number( data, b'V_hdr_hot_min', 2 )
	ssc.data_set_number( data, b'N_max_hdr_diams', 10 )
	ssc.data_set_number( data, b'L_rnr_pb', 25 )
	ssc.data_set_number( data, b'L_rnr_per_xpan', 70 )
	ssc.data_set_number( data, b'L_xpan_hdr', 20 )
	ssc.data_set_number( data, b'L_xpan_rnr', 20 )
	ssc.data_set_number( data, b'Min_rnr_xpans', 1 )
	ssc.data_set_number( data, b'northsouth_field_sep', 20 )
	ssc.data_set_number( data, b'N_hdr_per_xpan', 2 )
	ssc.data_set_number( data, b'offset_xpan_hdr', 1 )
	K_cpnt = [[ 0.90000000000000002,   0,   0.19,   0,   0.90000000000000002,   -1,   -1,   -1,   -1,   -1,   -1 ], [ 0,   0.59999999999999998,   0.050000000000000003,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.41999999999999998,   0,   0.14999999999999999 ], [ 0.050000000000000003,   0,   0.41999999999999998,   0,   0.59999999999999998,   0,   0.59999999999999998,   0,   0.14999999999999999,   0.59999999999999998,   0 ], [ 0.90000000000000002,   0,   0.19,   0,   0.90000000000000002,   -1,   -1,   -1,   -1,   -1,   -1 ]];
	ssc.data_set_matrix( data, b'K_cpnt', K_cpnt );
	D_cpnt = [[ 0.085000000000000006,   0.063500000000000001,   0.085000000000000006,   0.063500000000000001,   0.085000000000000006,   -1,   -1,   -1,   -1,   -1,   -1 ], [ 0.085000000000000006,   0.085000000000000006,   0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.063500000000000001,   0.085000000000000006,   0.085000000000000006,   0.085000000000000006 ], [ 0.085000000000000006,   0.063500000000000001,   0.085000000000000006,   0.063500000000000001,   0.085000000000000006,   -1,   -1,   -1,   -1,   -1,   -1 ]];
	ssc.data_set_matrix( data, b'D_cpnt', D_cpnt );
	L_cpnt = [[ 0,   0,   0,   0,   0,   -1,   -1,   -1,   -1,   -1,   -1 ], [ 0,   0,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   1,   0 ], [ 0,   1,   0,   1,   0,   0,   0,   1,   0,   0,   0 ], [ 0,   0,   0,   0,   0,   -1,   -1,   -1,   -1,   -1,   -1 ]];
	ssc.data_set_matrix( data, b'L_cpnt', L_cpnt );
	Type_cpnt = [[ 0,   1,   0,   1,   0,   -1,   -1,   -1,   -1,   -1,   -1 ], [ 1,   0,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   2,   0 ], [ 0,   2,   0,   2,   0,   1,   0,   2,   0,   0,   1 ], [ 0,   1,   0,   1,   0,   -1,   -1,   -1,   -1,   -1,   -1 ]];
	ssc.data_set_matrix( data, b'Type_cpnt', Type_cpnt );
	ssc.data_set_number( data, b'custom_sf_pipe_sizes', 0 )
	sf_rnr_diams = [[ -1 ]];
	ssc.data_set_matrix( data, b'sf_rnr_diams', sf_rnr_diams );
	sf_rnr_wallthicks = [[ -1 ]];
	ssc.data_set_matrix( data, b'sf_rnr_wallthicks', sf_rnr_wallthicks );
	sf_rnr_lengths = [[ -1 ]];
	ssc.data_set_matrix( data, b'sf_rnr_lengths', sf_rnr_lengths );
	sf_hdr_diams = [[ -1 ]];
	ssc.data_set_matrix( data, b'sf_hdr_diams', sf_hdr_diams );
	sf_hdr_wallthicks = [[ -1 ]];
	ssc.data_set_matrix( data, b'sf_hdr_wallthicks', sf_hdr_wallthicks );
	sf_hdr_lengths = [[ -1 ]];
	ssc.data_set_matrix( data, b'sf_hdr_lengths', sf_hdr_lengths );
	ssc.data_set_number( data, b'tanks_in_parallel', 1 )
	ssc.data_set_number( data, b'has_hot_tank_bypass', 0 )
	ssc.data_set_number( data, b'T_tank_hot_inlet_min', 400 )
	ssc.data_set_number( data, b'tes_pump_coef', 0.14999999999999999 )
	ssc.data_set_number( data, b'V_tes_des', 1.8500000000000001 )
	ssc.data_set_number( data, b'custom_tes_p_loss', 0 )
	k_tes_loss_coeffs = [[ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ]];
	ssc.data_set_matrix( data, b'k_tes_loss_coeffs', k_tes_loss_coeffs );
	ssc.data_set_number( data, b'custom_tes_pipe_sizes', 0 )
	tes_diams = [[ -1 ]];
	ssc.data_set_matrix( data, b'tes_diams', tes_diams );
	tes_wallthicks = [[ -1 ]];
	ssc.data_set_matrix( data, b'tes_wallthicks', tes_wallthicks );
	tes_lengths = [[ 0,   90,   100,   120,   0,   30,   90,   80,   80,   120,   80 ]];
	ssc.data_set_matrix( data, b'tes_lengths', tes_lengths );
	ssc.data_set_number( data, b'DP_SGS', 0 )
	ssc.data_set_number( data, b'adjust:constant', 4 )
	module = ssc.module_create(b'trough_physical')	
	ssc.module_exec_set_print( 0 );
	if ssc.module_exec(module, data) == 0:
		print ('trough_physical simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" );
	ssc.module_free(module)
	ssc.data_set_number( data, b'system_use_lifetime_output', 0 )
	ssc.data_set_number( data, b'analysis_period', 25 )
	ssc.data_set_number( data, b'enable_interconnection_limit', 0 )
	ssc.data_set_number( data, b'grid_interconnection_limit_kwac', 20000 )
	ssc.data_set_array_from_csv( data, b'grid_curtailment', b'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/PTcase/grid_curtailment.csv');
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
	ssc.data_set_number( data, b'ppa_soln_mode', 0 )
	ppa_price_input =[ 0.13 ];
	ssc.data_set_array( data, b'ppa_price_input',  ppa_price_input);
	ssc.data_set_number( data, b'ppa_escalation', 1 )
	federal_tax_rate =[ 21 ];
	ssc.data_set_array( data, b'federal_tax_rate',  federal_tax_rate);
	state_tax_rate =[ 7 ];
	ssc.data_set_array( data, b'state_tax_rate',  state_tax_rate);
	ssc.data_set_number( data, b'property_tax_rate', 0 )
	ssc.data_set_number( data, b'prop_tax_cost_assessed_percent', 100 )
	ssc.data_set_number( data, b'prop_tax_assessed_decline', 0 )
	ssc.data_set_number( data, b'real_discount_rate', 6.4000000000000004 )
	ssc.data_set_number( data, b'inflation_rate', 2.5 )
	ssc.data_set_number( data, b'insurance_rate', 0.5 )
	ssc.data_set_number( data, b'system_capacity', 99900 )
	om_fixed =[ 0 ];
	ssc.data_set_array( data, b'om_fixed',  om_fixed);
	ssc.data_set_number( data, b'om_fixed_escal', 0 )
	om_production =[ 4 ];
	ssc.data_set_array( data, b'om_production',  om_production);
	ssc.data_set_number( data, b'om_production_escal', 0 )
	om_capacity =[ 66 ];
	ssc.data_set_array( data, b'om_capacity',  om_capacity);
	ssc.data_set_number( data, b'om_capacity_escal', 0 )
	om_fuel_cost =[ 0 ];
	ssc.data_set_array( data, b'om_fuel_cost',  om_fuel_cost);
	ssc.data_set_number( data, b'om_fuel_cost_escal', 0 )
	om_replacement_cost1 =[ 0 ];
	ssc.data_set_array( data, b'om_replacement_cost1',  om_replacement_cost1);
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
	ptc_fed_amount =[ 0 ];
	ssc.data_set_array( data, b'ptc_fed_amount',  ptc_fed_amount);
	ssc.data_set_number( data, b'ptc_fed_term', 10 )
	ssc.data_set_number( data, b'ptc_fed_escal', 0 )
	ptc_sta_amount =[ 0 ];
	ssc.data_set_array( data, b'ptc_sta_amount',  ptc_sta_amount);
	ssc.data_set_number( data, b'ptc_sta_term', 10 )
	ssc.data_set_number( data, b'ptc_sta_escal', 0 )
	ssc.data_set_number( data, b'depr_alloc_macrs_5_percent', 90 )
	ssc.data_set_number( data, b'depr_alloc_macrs_15_percent', 1.5 )
	ssc.data_set_number( data, b'depr_alloc_sl_5_percent', 0 )
	ssc.data_set_number( data, b'depr_alloc_sl_15_percent', 2.5 )
	ssc.data_set_number( data, b'depr_alloc_sl_20_percent', 3 )
	ssc.data_set_number( data, b'depr_alloc_sl_39_percent', 0 )
	ssc.data_set_number( data, b'depr_alloc_custom_percent', 0 )
	depr_custom_schedule =[ 0 ];
	ssc.data_set_array( data, b'depr_custom_schedule',  depr_custom_schedule);
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
	degradation =[ 0 ];
	ssc.data_set_array( data, b'degradation',  degradation);
	roe_input =[ 0 ];
	ssc.data_set_array( data, b'roe_input',  roe_input);
	ssc.data_set_number( data, b'loan_moratorium', 0 )
	ssc.data_set_number( data, b'system_use_recapitalization', 0 )
	ssc.data_set_number( data, b'total_installed_cost', 562201536 )
	ssc.data_set_number( data, b'salvage_percentage', 0 )
	ssc.data_set_number( data, b'construction_financing_cost', 28110076 )
	ssc.data_set_number( data, b'depr_stabas_method', 1 )
	ssc.data_set_number( data, b'depr_fedbas_method', 1 )
	ssc.data_set_number( data, b'cp_capacity_payment_esc', 0 )
	ssc.data_set_number( data, b'cp_capacity_payment_type', 0 )
	cp_capacity_payment_amount =[ 0 ];
	ssc.data_set_array( data, b'cp_capacity_payment_amount',  cp_capacity_payment_amount);
	cp_capacity_credit_percent =[ 0 ];
	ssc.data_set_array( data, b'cp_capacity_credit_percent',  cp_capacity_credit_percent);
	ssc.data_set_number( data, b'cp_system_nameplate', 99.900001525878906 )
	ssc.data_set_number( data, b'cp_battery_nameplate', 0 )
	grid_curtailment_price =[ 0 ];
	ssc.data_set_array( data, b'grid_curtailment_price',  grid_curtailment_price);
	ssc.data_set_number( data, b'grid_curtailment_price_esc', 0 )
	module = ssc.module_create(b'singleowner')	
	ssc.module_exec_set_print( 0 );
	if ssc.module_exec(module, data) == 0:
		print ('singleowner simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" );
	ssc.module_free(module)
	print ('Weather File = ', weather)
	annual_energy = ssc.data_get_number(data, b'annual_energy');
	print ('Annual Net Electrical Energy Production = ', annual_energy)
	annual_thermal_consumption = ssc.data_get_number(data, b'annual_thermal_consumption');
	#print ('Annual Freeze Protection = ', annual_thermal_consumption)
	annual_tes_freeze_protection = ssc.data_get_number(data, b'annual_tes_freeze_protection');
	#print ('Annual TES Freeze Protection = ', annual_tes_freeze_protection)
	annual_field_freeze_protection = ssc.data_get_number(data, b'annual_field_freeze_protection');
	#print ('Annual Field Freeze Protection = ', annual_field_freeze_protection)
	capacity_factor = ssc.data_get_number(data, b'capacity_factor');
	#print ('Capacity factor = ', capacity_factor)
	annual_W_cycle_gross = ssc.data_get_number(data, b'annual_W_cycle_gross');
	#print ('Power cycle gross electrical output = ', annual_W_cycle_gross)
	kwh_per_kw = ssc.data_get_number(data, b'kwh_per_kw');
	#print ('First year kWh/kW = ', kwh_per_kw)
	conversion_factor = ssc.data_get_number(data, b'conversion_factor');
	#print ('Gross-to-net conversion = ', conversion_factor)
	annual_total_water_use = ssc.data_get_number(data, b'annual_total_water_use');
	#print ('Annual Water Usage = ', annual_total_water_use)
	ppa = ssc.data_get_number(data, b'ppa');
	#print ('PPA price (year 1) = ', ppa)
	lppa_nom = ssc.data_get_number(data, b'lppa_nom');
	#print ('Levelized PPA price (nominal) = ', lppa_nom)
	lppa_real = ssc.data_get_number(data, b'lppa_real');
	#Print ('Levelized PPA price (real) = ', lppa_real)
	lcoe_nom = ssc.data_get_number(data, b'lcoe_nom');
	#print ('Levelized COE (nominal) = ', lcoe_nom)
	lcoe_real = ssc.data_get_number(data, b'lcoe_real');
	#print ('Levelized COE (real) = ', lcoe_real)
	project_return_aftertax_npv = ssc.data_get_number(data, b'project_return_aftertax_npv');
	#print ('Net present value = ', project_return_aftertax_npv)
	flip_actual_irr = ssc.data_get_number(data, b'flip_actual_irr');
	#print ('Internal rate of return (IRR) = ', flip_actual_irr)
	flip_actual_year = ssc.data_get_number(data, b'flip_actual_year');
	#print ('Year IRR is achieved = ', flip_actual_year)
	project_return_aftertax_irr = ssc.data_get_number(data, b'project_return_aftertax_irr');
	#print ('IRR at end of project = ', project_return_aftertax_irr)
	cost_installed = ssc.data_get_number(data, b'cost_installed');
	#print ('Net capital cost = ', cost_installed)
	size_of_equity = ssc.data_get_number(data, b'size_of_equity');
	#print ('Equity = ', size_of_equity)
	size_of_debt = ssc.data_get_number(data, b'size_of_debt');
	#print ('Size of debt = ', size_of_debt)
	ssc.data_free(data);

	return([annual_energy, capacity_factor, ppa, lcoe_real, project_return_aftertax_npv, cost_installed])