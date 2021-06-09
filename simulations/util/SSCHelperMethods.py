#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:50:27 2021

@author: gabrielsoto
"""

import numpy as np
import pint
from pyomo.environ import units

class SSCHelperMethods(object):
    
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
        
        a = -1.0e-10*u.J/u.g/u.kelvin**4
        b = 2.0e-7*u.J/u.g/u.kelvin**3
        c = 5.0e-6*u.J/u.g/u.kelvin**2
        d = 1.4387*u.J/u.g/u.kelvin
        
        cp = (a*T**3 + b*T**2 + c*T + d).to('J/g/kelvin')
        
        return cp 


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