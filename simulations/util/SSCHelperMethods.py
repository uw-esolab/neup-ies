#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:50:27 2021

@author: gabrielsoto
"""

import numpy as np
import pint

class SSCHelperMethods(object):
    
    def define_unit_registry():
        u = pint.UnitRegistry()
        return u
        

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