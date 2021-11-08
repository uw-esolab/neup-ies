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


#parse the dispatch factors file
dispatch_data = []
with open ("../data/dispatch_factors_ts.csv") as f:
    for line in f:
        dispatch_data.append(float(line.strip()))
dispatch_data=np.array(dispatch_data)

# create empty dictionary for results
results={}

p_refs=[465,525,600,700,800]
cycle_max_fracs=[1.05]
remove_tes_costs=[False]
exaggerates=[1]#,1.5,2] 
for p_ref in p_refs:
    results[p_ref]={}
    
    for cycle_max_frac in cycle_max_fracs: #scoping indicates low sensitivity
        
        results[p_ref][cycle_max_frac]={}                    
        
        #try disabling the TES cost for some cases
        for remove_tes_cost in remove_tes_costs: #True
            
            #only do this option for first cycle_max_frac
            if cycle_max_frac != 1.05 and remove_tes_cost:
                continue
                
            results[p_ref][cycle_max_frac][remove_tes_cost]={}
            
            # make the PPA profile more extreme systematically
            for exaggerate in exaggerates:
                
                #load in the typical model1 simulation as starting point
                with open("../json/model1.json") as f:
                    base_sim = json.load(f)
                
                turbine_unit_cost = 225 #Cory email 8th Nov 2021
                turbine_ref = 950*base_sim["SSC_inputs"]["design_eff"]
                turbine_premium = turbine_unit_cost*(p_ref-turbine_ref)
                
                #add on the turbine premium - but then scale down to ensure that plant cost is linked to reactor rating not turbine
                base_sim["SSC_inputs"]["nuclear_spec_cost"]=(4500*turbine_ref+turbine_premium)/p_ref
                print(base_sim["SSC_inputs"]["nuclear_spec_cost"])
                
                if remove_tes_cost:
                    base_sim["SSC_inputs"]["tes_spec_cost"]=0.0
                else:
                    base_sim["SSC_inputs"]["tes_spec_cost"]=28.4 # lower with higher delta T
                    
                base_sim["SSC_inputs"]["cycle_max_frac"]=cycle_max_frac
                base_sim["SSC_inputs"]["P_ref"]=p_ref
                base_financing=base_sim["SSC_inputs"]["construction_financing_cost"]
                
                results[p_ref][cycle_max_frac][remove_tes_cost][exaggerate]=[]
                
                #increase the difference of the factors from unity
                new_dispatch_data = list((dispatch_data-1)*exaggerate+1)
                with open("../data/temp_dispatch_factors.csv","w") as f:
                    for item in new_dispatch_data:
                        f.write("{:.5f}\n".format(item))
                
                #sweep the TS Hours 
                for tshours in range(1): #6 if p_ref >465 else 1):
                    
                    base_sim["SSC_inputs"]["tshours"]=tshours
                    
                    #4 years, 7%. Built into initial financing estimate.
                    #Add on extra cost of big turbine and TES
                    yrs=4.0
                    rate=0.07
                    kWth=950000
                    extra_financing = yrs*rate*(1000*turbine_premium+kWth*tshours*base_sim["SSC_inputs"]["tes_spec_cost"])
                    base_sim["SSC_inputs"]["construction_financing_costs"]=base_financing+extra_financing
                    
                    #base_sim["SSC_inputs"]["construction_financing_costs"]=base_sim["SSC_inputs"]["total_installed_cost"]*yrs*rate 
                    
                    
                    with open("../json/tmp.json","w") as f:
                        json.dump(base_sim,f)
                    
                    # =============================================================================
                    #     Run Simulation
                    # =============================================================================
                    
                    # defining directories
                    nuctes = NuclearTES.NuclearTES(json_name="tmp",is_dispatch=True,log_dispatch_targets=False)
                    output_file = 'output_P{0:.0f}_C{1:.2f}_{2}_E{3:.1f}_T{4:.0f}.csv'.format(p_ref,cycle_max_frac,"T" if remove_tes_cost else "F",exaggerate,tshours)
                    print("attempting to run "+output_file)#Output results
                    nuctes.run_sim( run_loop=True, export=True, filename=output_file )
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
                    
                    results[p_ref][cycle_max_frac][remove_tes_cost][exaggerate].append(ppa)
                    
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

#Save results                           
with open("output_ppa_results.json","w") as f:
    json.dump(results,f)          

#Output results
for exaggerate in exaggerates:
    print("Exagerate "+str(exaggerate))
    for cycle_max_frac in cycle_max_fracs:
        print("Turbine Factor "+str(cycle_max_frac))
        print("Powers/TSHours,0,1,2,3,4,5")
        first_item=results[p_refs[0]][cycle_max_frac][False][exaggerate][0]
        for p_ref in p_refs:
            print("{0},{1}".format(str(p_ref),",".join(["{:.2f}".format(item/first_item) for item in results[p_ref][cycle_max_frac][False][exaggerate]])))

