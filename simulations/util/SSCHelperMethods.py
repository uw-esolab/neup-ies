#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:50:27 2021

@author: gabrielsoto
"""

import numpy as np
import pint
from pyomo.environ import units
from util.FileMethods import FileMethods
from scipy.optimize import curve_fit,fsolve
from iapws import IAPWS97

class SSCHelperMethods(object):
    """
    The SSCHelperMethods class is a util class for any auxilliary methods that 
    help us link up better with SSC or create better inputs for SSC. Typically
    for methods that can be used across all modules and dispatchers, so better
    to keep in one designated area. 
    """
    
    def define_unit_registry():
        """ Method to define unique unit registry from Pint
        
        This method creates a unique unit registry from the Pint Python package.
        This unit registry object is passed through classes to ensure that all
        use the same units.
        
        Outputs:
            u_pint (obj) : pint UnitRegistry object
            
        """
        # create unique unit registry
        u_pint = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
        
        # define currency units
        u_pint.define('USD = [currency]')
        u_pint.define('cents = 0.01 USD')
        
        # defining aliases
        u_pint.define('@alias USD = dollar')
        u_pint.define('@alias cents = cent')
        
        return u_pint
    

    def convert_to_pyomo_unit(params, param_str):
        """ Method to evaluate a Pyomo unit from a given Pint unit
        
        This method was created because Pyomo Environment units and the Pint 
        UnitRegistry use units in different ways and syntax. Pint seemed to have
        more functionality and use so I defaulted to using those in the generation
        and checking of Pyomo Parameters. Essentially what happens:
            (1) The entire `params` dictionary is taken as input as well as a single
               string for a specific parameter. 
            (2) A string representation of the unit(s) are created from the input
               e.g. u.m/u.s -> 'meter / second'
            (3) The algorithm below checks for any weirdness -e.g. are there multiple
               units, are they being multiplied/divided? etc.
            (4) The algorithm then converts the string representation of the Pint unit
               into the syntax for a Pyomo Environment unit.
            (5) If all goes well, it will return the evaluation of the Pyomo unit string,
               so on the Return end of the function call it will return a Pyomo unit. 
        
        Inputs:
            params (dict)    : dictionary of all Pyomo parameters to be used in optimization
            param_str (str)  : single string to be extracted from params dictionary
        Outputs:
            cmd (evaluated string)  : evaluation of a string called `unit_cmd` representing
                                      Pyomo unit syntax. If failed, return None.
            
        """
        
        u_pyomo = units # this line is necessary to carry out final evaluation
        Quant = params[param_str] # pint Quantity taken from params dictionary
    
        if hasattr(Quant,'u'):
            Unit  = Quant.u                     # Unit object of the previous pint Quantity
            unit_str   = Unit.format_babel()    # retrieving string representation of unit
            unit_cmd_start = 'u_pyomo.'         # start of command, calling pyomo unit package
    
            # empty list of unit strings
            clean_str_mult = []           # clean strings for multiplied units
            clean_str_div  = []           # clean strings for divided units
            str_contains_multipl  = False # flag for finding multiplication
            str_contains_division = False # flag for finding division during the multiplicative search
            unit_str_division = []        # empty, used to log division terms in multiplicative search
    
            # are any units multiplied? all multiplicative terms are automatically moved to the first terms
            if unit_str.find('*') != -1:
                
                # flag that multiplication exists
                str_contains_multipl = True
                
                # split unit str between multiplications
                unit_strs_mult = unit_str.split('*')
                
                # cycle through string split by '*'
                for m in range(len(unit_strs_mult)):
                    # for each split, making sure we don't grab blank spaces
                    clean_strs_mult = unit_strs_mult[m].split(' ')
                    # indeces of parts of string without blank spaces
                    clean_inds_mult  = [bool(len(clean_strs_mult[x]) != 0) for x in range(len(clean_strs_mult))]
                    # log necessary unit strings
                    clean_str_mult.append( np.extract(clean_inds_mult, clean_strs_mult)[0] )
                    
                    # checking if any of the strings contain division symbols
                    if unit_strs_mult[m].find('/') != -1:
                        str_contains_division = True            # flag for finding division during the multiplicative search
                        unit_str_division = unit_strs_mult[m]   # log division terms in multiplicative searchs
    
            # are any units divided? all divided terms are automatically moved to the last terms
            if unit_str.find('/') != -1:
                
                # use the previously found string with division in it, else just use the input
                unit_str_div = unit_str_division if str_contains_division else unit_str
                
                # flag that division exists
                str_contains_division = True 
                    
                # split unit str between multiplications
                unit_strs_div = unit_str_div.split('/')
                
                # cycle through string split by '/'
                for n in range(len(unit_strs_div)):
                    # for each split, making sure we don't grab blank spaces
                    clean_strs_div = unit_strs_div[n].split(' ')
                    # indeces of parts of string without blank spaces
                    clean_inds_div  = [bool(len(clean_strs_div[x]) != 0) for x in range(len(clean_strs_div))]
                    # grab extracted clean string
                    possible_str = np.extract(clean_inds_div,clean_strs_div)[0]
                    # log necessary unit strings
                    if possible_str in clean_str_mult:
                        # duplicative term
                        pass 
                    else:
                        # log term
                        clean_str_div.append( possible_str ) 
        
            # collecting all units
            unit_cmd = ''
            # multiplication
            if str_contains_multipl:
                unit_cmd  +=  unit_cmd_start + ('*'+unit_cmd_start).join(clean_str_mult)
            # division
            if str_contains_division:
                unit_cmd_div_start = '/' if str_contains_multipl else ''
                unit_cmd += unit_cmd_div_start+unit_cmd_start + ('/'+unit_cmd_start).join(clean_str_div)
            # neither
            if str_contains_division is False and str_contains_multipl is False:
                unit_cmd += unit_cmd_start + unit_str
            return eval(unit_cmd)
        else:
            return None
        

    def interpret_user_defined_cycle_data(ud_array):
        """ Method to return user defined cycle data (written by LORE team)
        
        Inputs:
            ud_array (list of list) : table of user defined data as nested lists
        Outputs:
            output_dict (dict) : dictionary of useful values from table
            
        """
        data = np.array(ud_array)
            
        i0 = 0
        nT = np.where(np.diff(data[i0::,0])<0)[0][0] + 1 
        Tpts = data[i0:i0+nT,0]
        mlevels = [data[j,1] for j in [i0,i0+nT,i0+2*nT]]
        
        i0 = 3*nT
        nm = np.where(np.diff(data[i0::,1])<0)[0][0] + 1 
        mpts = data[i0:i0+nm,1]
        Tamblevels = [data[j,2] for j in [i0,i0+nm,i0+2*nm]]
        
        i0 = 3*nT + 3*nm
        nTamb = np.where(np.diff(data[i0::,2])<0)[0][0] + 1 
        Tambpts = data[i0:i0+nTamb,2]
        Tlevels = [data[j,0] for j in [i0,i0+nm,i0+2*nm]]
        
        output_dict = {'nT':nT, 'Tpts':Tpts, 'Tlevels':Tlevels, 
                       'nm':nm, 'mpts':mpts, 'mlevels':mlevels, 
                       'nTamb':nTamb, 'Tambpts':Tambpts, 'Tamblevels':Tamblevels}
        
        return output_dict


    def get_cp_htf( u, T, rec_htf):
        """ Method to calculate specific heat of some heat transfer fluid (written by LORE team)
        
        Inputs:
            u (unitRegistry) : Pint unit registry
            T (float Quant)  : Temperature at which to find specific heat of HTF, in units of Kelvin
            rec_htf (int)    : integer value representing HTF in SSC table
        Outputs:
            cp (float Quant) : specific heat value in J/g/K
            
        """
        
        # HTF: 17 = Salt (60% NaNO3, 40% KNO3)
        if rec_htf != 17:
            print ('HTF %d not recognized'%rec_htf)
            return 0.0
        
        T = T.to('kelvin')
        
        a = -1.0e-10 * u.J/u.g/u.kelvin**4
        b =  2.0e-7  * u.J/u.g/u.kelvin**3
        c =  5.0e-6  * u.J/u.g/u.kelvin**2
        d =  1.4387  * u.J/u.g/u.kelvin
        
        cp = (a*T**3 + b*T**2 + c*T + d).to('J/g/kelvin')
        
        return cp 


    def get_rho_htf( u, T, rec_htf):
        """ Method to calculate density of some heat transfer fluid 
        
        This method calculates density at average temperature for a heat transfer
        fluid. In our case, we use molten salt. This method is copied from SSC,
        specifically the htf_props cpp method. 
        
        Inputs:
            u (unitRegistry) : Pint unit registry
            T (float Quant)  : Temperature at which to find specific heat of HTF, in units of Kelvin
            rec_htf (int)    : integer value representing HTF in SSC table
        Outputs:
            rho (float Quant) : fluid density value in kg/m^3
            
        """
        
        # HTF: 17 = Salt (60% NaNO3, 40% KNO3)
        if rec_htf != 17:
            print ('HTF %d not recognized'%rec_htf)
            return 0.0
        
        T = T.to('kelvin')

        a = -1.0e-07     * u.kg/u.m**3/u.kelvin**3
        b =  0.0002      * u.kg/u.m**3/u.kelvin**2
        c = -0.7875      * u.kg/u.m**3/u.kelvin
        d =  2299.4      * u.kg/u.m**3
        max_rho = 1000.0 *u.kg/u.m**3
        
        rho = (a*T**3 + b*T**2 + c*T + d).to('kg/m^3')
        rho = max(rho, max_rho)
        
        return rho 


    def get_visc_htf( u, T, rec_htf):
        """ Method to calculate viscosity of some heat transfer fluid (written by LORE team)
        
        Inputs:
            u (unitRegistry) : Pint unit registry
            T (float Quant)  : Temperature at which to find specific heat of HTF, in units of Kelvin
            rec_htf (int)    : integer value representing HTF in SSC table
        Outputs:
            visc (float Quant) : viscosity value in Ns/m^2
            
        """
        
        # HTF: 17 = Salt (60% NaNO3, 40% KNO3)
        if rec_htf != 17:
            print ('HTF %d not recognized'%rec_htf)
            return 0.0
        
        T = T.to('celsius')

        a = -1.473302e-10 * u.N*u.s/u.m**2/u.delta_degC**3
        b =  2.279989e-7  * u.N*u.s/u.m**2/u.delta_degC**2
        c = -1.199514e-4  * u.N*u.s/u.m**2/u.delta_degC
        d =  0.02270616   * u.N*u.s/u.m**2
        max_visc = 1e-4   * u.N*u.s/u.m**2
        
        visc = (a*T**3 + b*T**2 + c*T + d).to('N*s/m^2')
        visc = max(visc, max_visc)
        
        return visc 
    
    
    def get_ambient_T_corrections_from_udpc_inputs( u, Tamb, ud_array):
        """ Method to calculate specific heat of some heat transfer fluid (written by LORE team)
        
        Inputs:
            u (unitRegistry)         : Pint unit registry
            Tamb (float)             : ambient temperature, not units (intended to be C)
            ud_array (list of lists) : table of user defined data as nested lists
        Outputs:
            etamult (float) : multiplication correction factor for eta
            wmult (float)   : multiplication correction factor for w
            
        """
        # grabbing a dictionary of useful user-defined data for power cycle
        ud_dict = SSCHelperMethods.interpret_user_defined_cycle_data( ud_array )
        
        n = len(Tamb)  # Tamb = set of ambient temperature points for each dispatch time step
        
        Tambpts = np.array(ud_dict['Tambpts'])*u.degC
        i0 = 3*ud_dict['nT']+3*ud_dict['nm']+ud_dict['nTamb']  # first index in udpc data corresponding to performance at design point HTF T, and design point mass flow
        npts = ud_dict['nTamb']
        etapts = [ ud_array[j][3]/ud_array[j][4] for j in range(i0, i0+npts)]
        wpts = [ ud_array[j][5] for j in range(i0, i0+npts)]
        
        etamult  = np.ones(n)
        wmult = np.ones(n) 
        Tstep = Tambpts[1] - Tambpts[0]
        for j in range(n):
            i = max(0, min( int((Tamb[j] - Tambpts[0]) / Tstep), npts-2) )
            r = (Tamb[j] - Tambpts[i]) / Tstep
            etamult[j] = etapts[i] + (etapts[i+1] - etapts[i])*r
            wmult[j] = wpts[i] + (wpts[i+1] - wpts[i])*r

        return etamult, wmult


    def get_linearized_ud_params( ud_array, q_pb_design, SSC_dict):
        """ Method to calculate linearized user defined params (written by LORE team)
        
        Inputs:
            ud_array (list of lists)  : table of user defined data as nested lists
            q_pb_design (float Quant) : power block thermal power input rating, units of MW
            SSC_dict (dict)           : dictionary of SSC inputs
        Outputs:
            etap (float) : slope of linearized eta_p
            b (float)    : intercept of linearized eta_p
            
        """
        # grabbing a dictionary of useful user-defined data for power cycle
        ud_dict = SSCHelperMethods.interpret_user_defined_cycle_data( ud_array )
        
        # getting relevant eta points
        eta_adj_pts = [ud_array[p][3]/ud_array[p][4] for p in range(len(ud_array)) ]
        
        xpts = ud_dict['mpts']
        step = xpts[1] - xpts[0]
        
        # Interpolate for cycle performance at specified min/max load points
        fpts = [SSC_dict['cycle_cutoff_frac'], SSC_dict['cycle_max_frac']]
        q, eta = [ [] for v in range(2)]
        for j in range(2):
            p = max(0, min(int((fpts[j] - xpts[0]) / step), len(xpts)-2) )  # Find first point in user-defined array of load fractions for interpolation
            i = 3*ud_dict['nT'] + ud_dict['nm'] + p    # Index of point in full list of udpc points (at design point ambient T)
            eta_adj = eta_adj_pts[i] + (eta_adj_pts[i+1] - eta_adj_pts[i])/step * (fpts[j] - xpts[p])
            eta.append(eta_adj * SSC_dict['design_eff'])
            q.append(fpts[j]*q_pb_design)
        
        # calculating eta_P as linear slope using max and min power/thermal output/input
        etap = (q[1]*eta[1]-q[0]*eta[0])/(q[1]-q[0])
        
        # intercept of linearized eta_P
        b = q[1]*(eta[1] - etap)
        
        return etap, b


    def estimate_receiver_pumping_parasitic( u, Tavg, dm_rec_design, SSC_dict, nonheated_length = 0.2):
        """ Method to calculate pumping parasitic of receiver (written by LORE team)
        
        Inputs:
            u (unitRegistry)            : Pint unit registry
            Tavg (float Quant)          : average temperature of plant, K
            dm_rec_design (float Quant) : power block thermal power input rating, units of MW
            SSC_dict (dict)             : dictionary of SSC inputs
            nonheated_length (float)    : non-heated length of total receiver height
        Outputs:
            etap (float) : slope of linearized eta_p
            b (float)    : intercept of linearized eta_p
            
        """
        
        # flow parameters
        rho  = SSCHelperMethods.get_rho_htf(  u, Tavg, SSC_dict['rec_htf'] )
        visc = SSCHelperMethods.get_visc_htf( u, Tavg, SSC_dict['rec_htf'] )

        # some design parameters
        D_rec    = SSC_dict['D_rec'] * u.m
        N_panels = SSC_dict['N_panels']
        d_tube_out = SSC_dict['d_tube_out'] * u.mm
        th_tube    = SSC_dict['th_tube'] * u.mm
        rec_height = SSC_dict['rec_height'] * u.m
        h_tower    = SSC_dict['h_tower'] * u.m
        eta_pump   = SSC_dict['eta_pump'] 
        g0 = 9.8 * u.m / u.s**2
        
        # number of flow paths
        npath = 1
        nperpath = SSC_dict['N_panels']
        if SSC_dict['Flow_type'] == 1 or SSC_dict['Flow_type'] == 2:
            npath = 2
            nperpath = int(SSC_dict['N_panels']/2)
        elif SSC_dict['Flow_type'] == 9:
            npath = int(SSC_dict['N_panels']/2)
            nperpath = 2
        
        # ============================
        # fluid mechanics of the pump
        
        # tube dimensions
        ntube = int(np.pi * D_rec / N_panels / d_tube_out)  # Number of tubes per panel
        m_per_tube = dm_rec_design / npath / ntube  # kg/s per tube
        tube_id = d_tube_out - 2*th_tube            # Tube ID in m
        Ac = ( 0.25*np.pi*(tube_id**2) ).to('mm^2')
        
        # flow velocities and fluid mechanics
        vel = (m_per_tube / rho / Ac).to('mm/s')  # HTF velocity
        Re = (rho * vel * tube_id / visc).to('')  # Reynolds number
        EoverD = (4.6e-5*u.m) / tube_id     # roughness factor (E/D)
        
        # friction factor
        ff = (-1.737*np.log(0.269*EoverD - 2.185/Re*np.log(0.269*EoverD+14.5/Re)))**-2
        fd = 4*ff 
        
        # total height and frictional pressure drop
        Htot = rec_height * (1 + nonheated_length)
        dp = 0.5*fd*rho*(vel**2) * (Htot/tube_id + 4*30 + 2*16) * nperpath  # Frictional pressure drop (Pa) (straight tube, 90deg bends, 45def bends)
        dp = dp.to('N/m^2')
        
        # gravitational pressure drop
        dp += rho * g0 * h_tower  # Add pressure drop from pumping up the tower
        if nperpath%2 == 1:   
            dp += rho * g0 * Htot  
        
        # pumpimg parasitic
        wdot = dp * dm_rec_design / rho / eta_pump   # Pumping parasitic at design point reciever mass flow rate (MWe)
        wdot = wdot.to('MW')
        
        return wdot
    

    def get_pc_persist_and_off_logs( param_dict, plant, npts ):
        """ Method to log the amount of time Power Cycle has been ON and OFF
        
        This method uses SSC output data from the previous run to log how long
        the Power Cycle has been both ON and OFF. One of the two outputs will be
        populated in this method, and there are a bunch of logic statements to
        correctly log the respective length of time. Method adapted from LORE Team. 
        
        TODO: can we just input another dictionary instead of passing the full Plant?
        
        Inputs:
            param_dict (dict)    : dictionary of Pyomo dispatch parameters
            plant (obj)          : the full PySAM Plant object. 
            npts (int)           : length of the SSC horizon
        Outputs:
            disp_pc_persist0 (int) : length of time PC has been ON in the past segment
            disp_pc_off0 (int)     : length of time PC has been OFF in the past segment
        """
        
        # cycle state before start of most recent set of simulation calls
        previous_pc_state = plant.SystemControl.pc_op_mode_initial
        # cycle state after most recent set of simulation calls
        tmp_final_op_mode = plant.Outputs.pc_op_mode_final
        current_pc_state  = tmp_final_op_mode if type(tmp_final_op_mode) == float \
                                else tmp_final_op_mode[npts-1]
        # times when cycle is not generating power
        tslice = slice(0,npts,1)
        is_pc_not_on = np.array( plant.Outputs.P_cycle[tslice] ) <= 1.e-3
        
        ###=== Persist Log ===### 
        # if PC is ON
        if current_pc_state == 1:
            # array of times (PC was generating power == True)
            is_pc_current = np.array( plant.Outputs.P_cycle[tslice] ) > 1.e-3 
            
        # if PC is STANDBY
        elif current_pc_state == 2:
            # array of times (PC was generating power == False) + (PC getting input energy == True) + (PC using startup power == False)
            is_pc_current = np.logical_and( \
                                np.logical_and( \
                                    np.array( plant.Outputs.P_cycle[tslice] ) <= 1.e-3, np.array( plant.Outputs.q_pb[tslice] ) >= 1.e-3 ), \
                                    np.array( plant.Outputs.q_dot_pc_startup[tslice] ) <= 1.e-3 )
        
        # if PC is STARTUP
        elif current_pc_state == 0:
            # array of times (PC using startup power == True)
            is_pc_current = np.array( plant.Outputs.q_dot_pc_startup[tslice] ) > 1.e-3
        
        # if PC is OFF
        elif current_pc_state == 3:
            # array of times (PC getting input energy + PC using startup power == False)
            is_pc_current = (np.array( plant.Outputs.q_dot_pc_startup[tslice] ) + np.array( plant.Outputs.q_pb[tslice] ) ) <= 1.e-3
        
        ###=== Indexing ===###
        ssc_time_step = 1   # 1 hour per time step
        n = npts            # length of ssc horizon
        
        ###=== OFF Log ===###
        # if PC is ON
        if current_pc_state == 1: 
            # returning 0 for OFF log
            disp_pc_off0 = 0.0
        
        # if PC is OFF for full simulated horizon
        elif is_pc_not_on.sum() == n:  
            # add all OFF positions in this current horizon to existing OFF log
            disp_pc_off0 = param_dict['Yd0'].to('hr').m + n*ssc_time_step  
        
        # if PC is OFF for some portion of current horizon
        else:
            # find indeces of changed OFF state
            i = np.where(np.abs(np.diff(is_pc_not_on)) == 1)[0][-1]
            # use index to find length of times PC was oFF
            disp_pc_off0 = int(n-1-i)*ssc_time_step         
        
        ###=== Final Indexing and Logging ===###
        # Plant has not changed state over this simulation window:
        if n == 1 or np.abs(np.diff(is_pc_current)).max() == 0:  
            # adding to existing persist array from Dispatch Params dictionary if state continued
            disp_pc_persist0 = n*ssc_time_step if previous_pc_state != current_pc_state else param_dict['Yu0'].to('hr').m + n*ssc_time_step
        # Plant *has* changed state over this simulation window:
        else:
            # find indeces of changed state
            i = np.where(np.abs(np.diff(is_pc_current)) == 1)[0][-1]
            # use index to find length of times PC was ON
            disp_pc_persist0 = int(n-1-i)*ssc_time_step
        
        return disp_pc_persist0, disp_pc_off0


    def linearize_indirectTES_eff(propn_LFR_flow = 0.8):
        raw_data = \
        """TF Temp	HTF mdot	Ambient T	Wdot cycle	Heat In
        500	0.300	30	0.237641861	0.31466289
        504.2	0.300	30	0.242543393	0.318919287
        508.4	0.300	30	0.247479625	0.323186151
        512.6	0.300	30	0.252450559	0.327456504
        516.8	0.300	30	0.257447519	0.331737323
        521.1	0.300	30	0.262470504	0.336021631
        525.3	0.300	30	0.267528191	0.340312917
        529.5	0.300	30	0.272611903	0.34461118
        533.7	0.300	30	0.277721641	0.348916421
        537.9	0.300	30	0.28286608	0.35322864
        542.1	0.300	30	0.288036545	0.357544348
        546.3	0.300	30	0.293241711	0.361867033
        550.5	0.300	30	0.298464228	0.366196696
        554.7	0.300	30	0.303712771	0.370533337
        558.9	0.300	30	0.308996014	0.374873467
        563.2	0.300	30	0.314305284	0.379220574
        567.4	0.300	30	0.319631904	0.38357117
        571.6	0.300	30	0.3250019	0.387932233
        575.8	0.300	30	0.330389247	0.392293296
        580	0.300	30	0.335811295	0.396664825
        500	1.000	30	0.766426017	0.810896011
        504.2	1.000	30	0.782518656	0.824028044
        508.4	1.000	30	0.798672022	0.837187987
        512.6	1.000	30	0.81489479	0.850382818
        516.8	1.000	30	0.83118696	0.863602072
        521.1	1.000	30	0.847548532	0.876849236
        525.3	1.000	30	0.863962156	0.890120822
        529.5	1.000	30	0.880427832	0.903420319
        533.7	1.000	30	0.896954234	0.916744238
        537.9	1.000	30	0.913524013	0.930092579
        542.1	1.000	30	0.930137169	0.943461853
        546.3	1.000	30	0.946802376	0.956855549
        550.5	1.000	30	0.963493608	0.970266689
        554.7	1.000	30	0.980228218	0.98370574
        558.9	1.000	30	0.996980177	0.997158747
        563.2	1.000	30	1.013766838	1.010632686
        567.4	1.000	30	1.030562174	1.02412407
        571.6	1.000	30	1.047383536	1.037636387
        575.8	1.000	30	1.064213574	1.051162659
        580	1.000	30	1.081052286	1.064706376
        500	1.200	30	0.889051057	0.924098734
        504.2	1.200	30	0.907668201	0.939515964
        508.4	1.200	30	0.926337397	0.954957615
        512.6	1.200	30	0.94504997	0.970434154
        516.8	1.200	30	0.963788568	0.985938605
        521.1	1.200	30	0.982561867	1.001477943
        525.3	1.200	30	1.001361192	1.017034727
        529.5	1.200	30	1.020169193	1.032622909
        533.7	1.200	30	1.038985868	1.048232025
        537.9	1.200	30	1.057802544	1.063865563
        542.1	1.200	30	1.076593194	1.079513056
        546.3	1.200	30	1.095375169	1.095184971
        550.5	1.200	30	1.114122442	1.110867353
        554.7	1.200	30	1.132835014	1.126570668
        558.9	1.200	30	1.15148686	1.14228096
        563.2	1.200	30	1.170077978	1.158008697
        567.4	1.200	30	1.188591019	1.173739923
        571.6	1.200	30	1.207025983	1.189481615
        575.8	1.200	30	1.225348168	1.205223307
        580	1.200	30	1.243566249	1.220971977
        560	0.300	0	0.338665461	0.375968966
        560	0.347	0	0.396304002	0.425821147
        560	0.395	0	0.452910184	0.474173122
        560	0.442	0	0.508414606	0.521171423
        560	0.490	0	0.562947396	0.567074225
        560	0.537	0	0.616161545	0.611700108
        560	0.584	0	0.668221882	0.655261892
        560	0.632	0	0.719102381	0.697829354
        560	0.679	0	0.768976549	0.739580424
        560	0.726	0	0.817532074	0.780323217
        560	0.774	0	0.864916436	0.820221708
        560	0.821	0	0.910426939	0.859314275
        560	0.868	0	0.9546188	0.89774396
        560	0.916	0	0.997622148	0.935311899
        560	0.963	0	1.039532411	0.97214718
        560	1.010	0	1.080297537	1.008274224
        560	1.058	0	1.119978254	1.043808163
        560	1.105	0	1.158279602	1.078543156
        560	1.153	0	1.195236281	1.112590845
        560	1.200	0	1.230761541	1.14595472
        560	0.300	30	0.310323332	0.375958499
        560	0.347	30	0.363129742	0.425817658
        560	0.395	30	0.415155378	0.474176611
        560	0.442	30	0.466356863	0.521185378
        560	0.490	30	0.51685565	0.567095158
        560	0.537	30	0.566348105	0.611728019
        560	0.584	30	0.614947006	0.655303758
        560	0.632	30	0.66266103	0.697885175
        560	0.679	30	0.709602953	0.73965369
        560	0.726	30	0.755469142	0.780410438
        560	0.774	30	0.800398402	0.820326374
        560	0.821	30	0.844347355	0.859436385
        560	0.868	30	0.887420105	0.897883514
        560	0.916	30	0.929330368	0.935465408
        560	0.963	30	0.970156221	0.972318133
        560	1.010	30	1.009854288	1.008462622
        560	1.058	30	1.048476621	1.044010516
        560	1.105	30	1.085728261	1.078762953
        560	1.153	30	1.121643907	1.112831576
        560	1.200	30	1.156145484	1.146212895
        560	0.300	55	0.233850765	0.375951522
        560	0.347	55	0.276411674	0.42581068
        560	0.395	55	0.318894505	0.474169633
        560	0.442	55	0.361195156	0.52118189
        560	0.490	55	0.403330976	0.567095158
        560	0.537	55	0.444989656	0.611734997
        560	0.584	55	0.486205897	0.655317714
        560	0.632	55	0.526944997	0.697909597
        560	0.679	55	0.567250334	0.739688579
        560	0.726	55	0.606826947	0.78046626
        560	0.774	55	0.64576159	0.820403128
        560	0.821	55	0.683993536	0.859534072
        560	0.868	55	0.721748341	0.898009113
        560	0.916	55	0.758652969	0.935618918
        560	0.963	55	0.794672719	0.972499554
        560	1.010	55	0.829738188	1.008671953
        560	1.058	55	0.863858053	1.044254736
        560	1.105	55	0.896763378	1.079035084
        560	1.153	55	0.928471515	1.113135106
        560	1.200	55	0.958869686	1.146551313
        500	1.000	0	0.830649092	0.810742502
        500	1.000	2.895	0.830649092	0.810742502
        500	1.000	5.789	0.830649092	0.810742502
        500	1.000	8.684	0.830649092	0.810742502
        500	1.000	11.58	0.830067849	0.810745991
        500	1.000	14.47	0.827508642	0.810759946
        500	1.000	17.37	0.822468306	0.810780879
        500	1.000	20.26	0.814435	0.810805301
        500	1.000	23.16	0.803200516	0.810833212
        500	1.000	26.05	0.789146568	0.810857634
        500	1.000	28.95	0.77277632	0.810885545
        500	1.000	31.84	0.7549226	0.810913456
        500	1.000	34.74	0.736123275	0.810934389
        500	1.000	37.63	0.717089717	0.810955322
        500	1.000	40.53	0.698108211	0.810976255
        500	1.000	43.42	0.67959517	0.81099021
        500	1.000	46.32	0.661628672	0.811000677
        500	1.000	49.21	0.644434273	0.811004166
        500	1.000	52.11	0.627942571	0.811004166
        500	1.000	55	0.612283696	0.810993699
        560	1.000	0	1.071379353	1.000344067
        560	1.000	2.895	1.071379353	1.000344067
        560	1.000	5.789	1.071379353	1.000344067
        560	1.000	8.684	1.071379353	1.000344067
        560	1.000	11.58	1.070589902	1.000347556
        560	1.000	14.47	1.067874541	1.000365
        560	1.000	17.37	1.062539245	1.000385933
        560	1.000	20.26	1.053968072	1.000413844
        560	1.000	23.16	1.041831359	1.000445244
        560	1.000	26.05	1.026432743	1.000480132
        560	1.000	28.95	1.008275388	1.000515021
        560	1.000	31.84	0.988226823	1.00054642
        560	1.000	34.74	0.966902992	1.000581309
        560	1.000	37.63	0.945084671	1.000612709
        560	1.000	40.53	0.92313622	1.000640619
        560	1.000	43.42	0.901578157	1.00066853
        560	1.000	46.32	0.880514584	1.000692952
        560	1.000	49.21	0.860240462	1.000710396
        560	1.000	52.11	0.840703738	1.000724352
        560	1.000	55	0.822069244	1.000731329
        580	1.000	0	1.153239266	1.064511
        580	1.000	2.895	1.153239266	1.064511
        580	1.000	5.789	1.153239266	1.064511
        580	1.000	8.684	1.153239266	1.064511
        580	1.000	11.58	1.152406439	1.064517978
        580	1.000	14.47	1.149621675	1.064531933
        580	1.000	17.37	1.144173601	1.064556355
        580	1.000	20.26	1.135402896	1.064587755
        580	1.000	23.16	1.122953874	1.064619155
        580	1.000	26.05	1.107112818	1.064657532
        580	1.000	28.95	1.088391571	1.06469242
        580	1.000	31.84	1.067648983	1.064730798
        580	1.000	34.74	1.045527027	1.064765686
        580	1.000	37.63	1.022841178	1.064800575
        580	1.000	40.53	0.999964473	1.064835463
        580	1.000	43.42	0.977443454	1.064866863
        580	1.000	46.32	0.9553909	1.064894774
        580	1.000	49.21	0.934119121	1.064915707
        580	1.000	52.11	0.913602091	1.06493664
        580	1.000	55	0.89400464	1.064947107"""
        
        #This function is used to parse the above table
        def get_x_and_PQ(raw_data):
            lines = raw_data.split("\n")[1:]
            
            y1=[] #W
            y2=[] #Q
            x=[[],[],[]] #Tf,mdot,Tamb
            
            for line in lines:
                items = line.strip().split()
                for j in range(3):
                    x[j].append(float(items[j]))
                y1.append(float(items[3]))
                y2.append(float(items[4]))
            return x,y1,y2
        
        #parse the Hamilton data
        x,raw_P,raw_Q= get_x_and_PQ(raw_data)
        
        #Turn Hamilton Data into Numpy arrays          
        x=np.array(x)
        raw_P=np.array(raw_P)
        raw_Q=np.array(raw_Q)
        
        #Add 73K to the hot temperature to match the Westinghouse design point of 633C
        x[0]+=73 
        
        def fn(x,a,b,c,d,
               e,
               f,
               g,
               h,
               i,
               j
               ):
            
            result=  a+b*x[0]+c*x[1]+d*x[2]
            result+= e*x[0]*x[1]
            result+= f*x[1]*x[2]
            result+= g*x[0]*x[2]
            result+= h*x[0]**2
            result+= i*x[1]**2
            result+=j*x[2]**2
            return result
        
        def find_popt(x,y,fn):
            popt,pcov=curve_fit(fn,x,y)
            predict=fn(x,*popt)
            residual=np.sum((y-predict)**2)
            # print(residual)
            return predict,popt
        
        predict_P,popt_P = find_popt(x,raw_P,fn)
        predict_Q,popt_Q = find_popt(x,raw_Q,fn)
        
        P=[]
        Q=[]
        for j,T_CSP in enumerate([550,633]): 
        
            #Proportion of turbine rating delivered by LFR only at full power . Can change   
            T_amb = 32.44 #C Ambient
            wec_P = 465.0 #MW LFR electrical power
            design_P = wec_P/propn_LFR_flow #MW 
            T_LFR = 633 #C source: Westinghouse PEPSE
            pressure = 33 #MPa source: westinghouse
            
            mdot_lfr_only = wec_P/design_P #mdot at which Salt2Steam kicks in
            
            #Steam properties
            lfr_steam = IAPWS97(T=T_LFR+273.15,P=pressure)
            csp_steam = IAPWS97(T=T_CSP+273.15,P=pressure)
            
            y=[[],[],[]] #Tf,mdot,T_amb
            
            #Vary mdot to get different P and Q
            mdots = [0.01]+[0.1*j for j in range(1,20)]
            for mdot in mdots:
                #calculate the temperature of the hot stream
                y[1].append(mdot)
                y[2].append(T_amb)
                if mdot<=mdot_lfr_only:
                    y[0].append(T_LFR)
                else:
                    #Mass-flow weighted enthalpy to turbine - then get temperature
                    h = (lfr_steam.h*mdot_lfr_only+csp_steam.h*(mdot-mdot_lfr_only))/mdot
                    T_mix = IAPWS97(h=h,P=pressure).T-273.15
                    y[0].append(T_mix)
            
        
            y = np.array(y)
            SAM_gen1 = fn(y,*popt_P)
            SAM_gen2 = fn(y,*popt_Q)
            Q.append(SAM_gen2)
            P.append(SAM_gen1)
        
        Q=np.array(Q)
        P=np.array(P) 
        
        def plot_P_and_Q(P,Q):

            #Plot a straight line for each case going through mdot=100%
            m = np.divide(P[:,11]-P[:,9],Q[:,11]-Q[:,9])
            c = P[:,10]-np.multiply(m,Q[:,10])
        
            return m,c
        
        m,c=plot_P_and_Q(P,Q)
        
        Q2=np.array(Q)
        lfr_steam = IAPWS97(T=633+273.15,P=pressure)
        lfr_steam_s = lfr_steam.s
        
        def turbine_outlet(p):
            #find pressure ratio that gives turbine isentropic efficiency of 90% and outlet T=570C
            p=p[0]
            s2s_steam=IAPWS97(T=570+273.15,P=p)
            delta_h_s=(lfr_steam.h-s2s_steam.h)/0.9
            s2s_steam_s=IAPWS97(h=lfr_steam.h-delta_h_s,P=p).s
            return s2s_steam_s-lfr_steam_s
        
        p=fsolve(turbine_outlet,[24])[0]
        
        s2s_steam = IAPWS97(T=570+273.15,P=p) #enthalpy of steam going to salt2steam
        cond_steam = IAPWS97(T=306,P=0.5) #enthalpy of steam at turbine exhaust
        
        extra_enthalpy = (lfr_steam.h-s2s_steam.h)/(lfr_steam.h-cond_steam.h) #the enthalpy that was extracted upstream of steam2salt relative to design point enthalpy
        
        #Now add this to Q2
        for j in range(int(10*propn_LFR_flow)):
            Q2[:,j]+=(1-x[1][j])*extra_enthalpy #mass flow weighting (1-x[1][j] is the mass flow through the steam2salt HX)
        # print(lfr_steam.h,s2s_steam.h,cond_steam.h,extra_enthalpy)
        
        P2=np.array(P)
        for j in range(int(10*propn_LFR_flow)):
            P2[:,j]+=(1-y[1][j])*extra_enthalpy
        
        m2,c2=plot_P_and_Q(P2,Q2)
        
        return m2,c2
