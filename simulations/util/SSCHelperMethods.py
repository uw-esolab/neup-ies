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
        
        # create unique unit registry
        u_pint = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
        
        # define currency units
        u_pint.define('USD = [currency]')
        u_pint.define('cents = 0.01 USD')
        
        # defining aliases
        u_pint.define('@alias USD = dollar')
        u_pint.define('@alias cents = cent')
        
        return u_pint
    

    def get_pyomo_unit(params, param_str):
        
        u_pyomo = units
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
                unit_cmd += '/'+unit_cmd_start + ('/'+unit_cmd_start).join(clean_str_div)
            # neither
            if str_contains_division is False and str_contains_multipl is False:
                unit_cmd += unit_cmd_start + unit_str
            return eval(unit_cmd)
        else:
            return None
        

    def interpret_user_defined_cycle_data(ud_ind_od):
        """
        Method written by NREL Lore team
        """
        data = np.array(ud_ind_od)
            
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