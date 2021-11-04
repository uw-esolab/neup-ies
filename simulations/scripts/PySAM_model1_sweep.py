#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:47:58 2020

Most recently tested against PySAM 2.1.4

@author: frohro
"""

import os,sys
sys.path.append('..')
import modules.NuclearTES as NuclearTES
import matplotlib.pyplot as plt
import pint
import numpy as np
import pyomo.environ as pe
from pylab import rc
import json

rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

pid = os.getpid()
print("PID = ", pid)

#load in the typical model1 simulation as starting point
with open("../json/model1.json") as f:
    base_sim = json.load(f)

# create empty dictionary for results
results={}

#two turbine max fractions
for cycle_max_frac in [1.05,1.5]: 
    base_sim["SSC_inputs"]["cycle_max_frac"]=cycle_max_frac
    
    results[cycle_max_frac]={}
    
    #try disabling the TES cost for some cases
    for remove_tes_cost in [False,True]:
        
        #only do this option for first cycle_max_frac
        if cycle_max_frac != 1.05 and remove_tes_cost:
            continue
        
        if remove_tes_cost:
            base_sim["SSC_inputs"]["tes_spec_cost"]=0.0
            
        results[cycle_max_frac][remove_tes_cost]={}
        
        # make the PPA profile more extreme systematically
        for exaggerate in [1,1.5,2.0]:
            results[cycle_max_frac][remove_tes_cost][exaggerate]=[]
            
            #increase the difference of the factors from unity
            for j in range(1,9):
                this_fac = "dispatch_factor"+str(int(j))
                base_sim["SSC_inputs"][this_fac]=(base_sim["SSC_inputs"][this_fac]-1)*exaggerate+base_sim["SSC_inputs"][this_fac]
            
            #sweep the TS Hours (nb currently no difference beyond 5)
            for tshours in [0,5,10,15]:
                if tshours == 0 and remove_tes_cost:
                    # duplicate case
                    results[cycle_max_frac][True][exaggerate].append(results[False][exaggerate][0])
                    continue
                
                base_sim["SSC_inputs"]["tshours"]=tshours
                
                with open("../json/tmp.json","w") as f:
                    json.dump(base_sim,f)
                
                # =============================================================================
                #     Run Simulation
                # =============================================================================
                
                # defining directories
                nuctes = NuclearTES.NuclearTES(json_name="tmp",is_dispatch=True,log_dispatch_targets=False)
                output_file = 'output.csv'
                nuctes.run_sim( run_loop=True, export=False, filename=output_file )
                nt = nuctes.Plant
                so = nuctes.SO
                
                print('Made it past execute.')
                
                # =============================================================================
                #     Display Results
                # =============================================================================
                
                annual_energy          = (nt.Outputs.annual_energy*u.kWh).to('TWh')
                capacity_factor        = nuctes.capacity_factor.magnitude * 100
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
                
                results[cycle_max_frac][remove_tes_cost][exaggerate].append(ppa)
                
                print('')
                print('                        Nuclear')
                print ('Annual energy (year 1)        =   ', annual_energy)
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
                
                # =============================================================================
                #     Plotting
                # =============================================================================
                
                from util.PostProcessing import Plots
                upl = Plots(nuctes, legend_offset = True, x_shrink=0.7)
                
                # 48 hour plot
                fig = plt.figure(figsize=[14,6])
                ax1 = fig.add_subplot(311)
                ax2 = fig.add_subplot(312)
                ax3 = fig.add_subplot(313)
                
                plt_allTime = False
                title = 'SSC Results - 48 hrs'
                start = 0
                end   = start + 72*2
            
                upl.plot_SSC_power_and_energy(ax1 , plot_all_time=plt_allTime, title_label=title, start_hr=start, end_hr=end, hide_x=True, x_legend=1.2, y_legend_L=1.0, y_legend_R=0.3)
                upl.plot_SSC_op_modes(ax2, plot_all_time=plt_allTime, start_hr=start, end_hr=end, hide_x=True )
                upl.plot_SSC_massflow(ax3, plot_all_time=plt_allTime, start_hr=start, end_hr=end, y_legend_L=0.8, y_legend_R=0.3)
                        
with open("output_ppa_results.json","w") as f:
    json.dump(results,f)          

