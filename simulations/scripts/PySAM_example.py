#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:47:58 2020

Most recently tested against PySAM 2.1.4

@author: frohro
"""

import modules.GenericSSCModule as GenericSSCModule

# defining directories
genmod = GenericSSCModule.GenericSSCModule()
genmod.run_sim()
nt = genmod.Plant
so = genmod.SO

print('Made it past execute.')
#print(gs.Outputs.export())  # as dictionary
