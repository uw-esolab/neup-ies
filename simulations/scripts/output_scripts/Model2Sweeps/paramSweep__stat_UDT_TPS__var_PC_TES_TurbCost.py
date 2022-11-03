#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 10:23:16 2022

@author: gabrielsoto
"""

import modules.DualPlantTES as DualPlantTES
from util.DualPlots import DualOutputExtraction
from util.FileMethods import FileMethods
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import numpy as np
import hashlib
from pylab import rc
import json as js
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# Static  - User Defined Tables + Tariff Price Schedule
# Varying - PC size, TES size
# =============================================================================

def get_turbine_cost( p ):
    """ Power law approximation to TOTAL cost
    """
    return 1630.418*p**(1-0.23)


# =============================================================================
# Parameters for Sweep
# =============================================================================
# possible different CSP design qdot. CURRENTLY HARDCODED SAM DEFAULT

sscH = 24  # 12 # 24
pyoH = 48  # 24 # 48
#json = 'model2_Hamilton_560_tariffx1_mod' # model2_Hamilton_560_tariffx1_mod model2_CAISO_Hamilton_mod

dispatch = True # True # False
run_loop = True

# costs and other params
nuc_spec_cost  = 4150-400 #remove the turbine cost from the nuclear cost (400 is from Cory)
tes_spec_cost  = 29.8
fin_yrs  = 4.0
fin_rate = 0.07



#TODO: first run the zero cases. Then redo the syn cases before procedding to opt as if these dont check out no point
model1=True
for json in ['model2_CAISO_Hamilton_mod']:#,'model2_Hamilton_560_tariffx1_mod']:
    
    for root_case in ['sweep']:#,'small1','small2','micro']:
        
        if root_case == 'zero':
            syn = False
        else:
            syn = True
        
        if syn:
            case='syn_'+root_case
        else:
            case=root_case
        
        if syn:
            solar_json="115"
            if model1:
                solar_json="1"
            if "CAISO" not in json:
                if root_case=='sweep':
                    q_dot_nuclear_des=950
                    p_cycle=np.array([655])
                    tshours=np.array([1.25])
                
                if root_case=='large':
                    q_dot_nuclear_des=1900
                    p_cycle=np.array([1105])
                    tshours=np.array([0.74])
    
                if root_case=='small1':
                    q_dot_nuclear_des=250
                    p_cycle=np.array([323])
                    tshours=np.array([2.54])
                    
                if root_case=='small2':
                    q_dot_nuclear_des=100
                    p_cycle=np.array([252])
                    tshours=np.array([3.25])
                    
                if root_case=='micro':
                    q_dot_nuclear_des=20
                    p_cycle=np.array([215])
                    tshours=np.array([3.82])
            else:
                if root_case=='sweep':
                    q_dot_nuclear_des=950
                    p_cycle=np.array([1060])
                    tshours=np.array([4.91])
                    if model1:
                        p_cycle=np.array([900])
                        tshours=np.array([4])
                
                if root_case=='large':
                    q_dot_nuclear_des=1900
                    p_cycle=np.array([1960])
                    tshours=np.array([4.49])
                    if model1:
                        p_cycle=np.array([1800])
                        tshours=np.array([4])
    
                if root_case=='small1':
                    q_dot_nuclear_des=250
                    p_cycle=np.array([383.7])
                    tshours=np.array([6.50])
                    
                if root_case=='small2':
                    q_dot_nuclear_des=100
                    p_cycle=np.array([249.5])
                    tshours=np.array([7.85])
                    
                if root_case=='micro':
                    q_dot_nuclear_des=20
                    p_cycle=np.array([174.7])
                    tshours=np.array([9.49])

        if case == "nuc":
            #nuclear only as solar should constantly defocus. Currently  hangs
            q_dot_nuclear_des=950
            solar_json = "115"
            p_cycle    = np.array([ 465]) 
            tshours    = np.array([ 0 ])
        
        elif case == "sol":
            #test case. should be consistent with roughyl 40% capacity factor for solar
            #at 15 c/kWh and nuclear at 450MWsh 6.7 c/kWh. Here 76% capacity factor gives
            #480 MWe produced. Estimated ppa = (0.4*115*15+450*6.5)/480=7.5 cf 7.7 (good enough)
            q_dot_nuclear_des=950
            solar_json = "115"
            p_cycle = np.array([465+115])
            tshours    = np.array([ 0 ])
        
        elif case == "sweep":
            q_dot_nuclear_des=950
            solar_json = "115"
            tshours    = np.array([ 0, 1,2,3,4,5,6,7,8,9,10 ])
            p_cycle    = np.array([1000,950,900,850,800,750,700,650,600,550]) 
            
        elif case == "small1":
            q_dot_nuclear_des=250
            solar_json = "115"
            tshours    = np.array([ 0,2,4,6,8,10,12,14,16 ])
            p_cycle    = np.array([600,525,450,375,300,225,150]) 
        
        elif case=='bigtest':
            #very large nuclear unit has been compared to Model1 
            #seems to produce sensible capacity factor and costs for default
            #Model1 seems to break down on capacity factor at extreme turbines but when
            #we set BOP equal we get within 0.1
            q_dot_nuclear_des=950e3
            solar_json="115"
            tshours=np.array([0])
            p_cycle=np.array([450e3])
            
        
        elif case=="large":
            q_dot_nuclear_des=1900
            solar_json = "115"
            tshours    = np.array([0,1,2,3,4,5,6,7,8,9,10])
            p_cycle    = np.array([1900,1800,1700,1600,1500,1400,1300,1200,1100,1000]) 
            
        elif case=="small2":
            q_dot_nuclear_des=100
            solar_json = "115"
            tshours    = np.array([ 0,2,4,6,8,10,12,14,16 ])
            p_cycle    = np.array([350,325,300,275,250,225,200,175,150,125]) 
        
        elif case=="micro":
            q_dot_nuclear_des=20
            solar_json="115"
            solar_json = "115"
            tshours    = np.array([ 0,2,4,6,8,10,12,14,16 ])
            p_cycle    = np.array([300,280,260,240,220,200,180,160,140,120]) 
        
        
        
        elif case=='zero':
            q_dot_nuclear_des=0.1
            solar_json = "115"
            tshours    = np.array([ 0,2,4,6,8,10,12,14,16 ])
            p_cycle    = np.array([220,205,190,175,160,145,130,115]) 
            
        elif case=="samtest":
            #gives 8.5 (default tariff) cf 8.6 in actual SAM for a specific test case (care on BOP cost). Good enough due to Gurobi/code differences
            #evidences correct turbine costs
            #have confirmed that nuclear specific costs go cf turbine by moving BOP costs here too and no difference in output
            q_dot_nuclear_des=0.1
            solar_json = "115"
            tshours    = np.array([2])
            p_cycle    = np.array([200]) 
        
    
    
        # TES sizes to sweep through
        
        # tshours    = np.array([ 7, 8, 9 ])
        # tshours    = np.array([ 0, 2, 4, 6, 8, 10 ])
        # tshours    = np.array([ 0, 2, 4, 6, 8, 10, 12, 14 ])
        # tshours    = np.array([ 0, 4, 8, 12, 16, 20 ])
        # tshours    = np.array([ 12, 14, 16, 18, 20, 22 ])
        
        # PC sizes to sweep through
        
        # p_cycle    = np.array([ 700,  650,  600]) 
        # p_cycle    = np.array([ 900, 850,  800,  750]) 
        # p_cycle    = np.array([ 1050, 1000,  950]) 
        
        # p_cycle    = np.array([ 750,  700,  650,  600 ]) 
        # p_cycle    = np.array([ 950,  900,  850,  800 ]) 
        # p_cycle    = np.array([ 1150, 1100, 1050, 1000 ]) 
        
        
        #load in solar parameters
        with open("../../../json/sol/"+solar_json+".json") as f:
            sol = js.load(f)
        
        
        #parameters which vary with system size
        sol_params=["focus_type","helio_area_tot","rec_height","D_rec","h_tower",
        "helio_positions","land_area_base","system_capacity","cp_system_nameplate"]
        #note:sol construction financing added later
        
        # =============================================================================
        # Initialize empty arrays
        # =============================================================================
        
        # formally defining the iterator arrays
        iterator1 = tshours  
        iterator2 = p_cycle
        
        # initializing an empty array
        empty   = np.zeros([  len(iterator1) , len(iterator2) ])
        empty2D = np.zeros([  len(iterator1) , len(iterator2), 2 ])
        
        # initializing output arrays
        defocus   = empty2D.copy()
        gen       = empty2D.copy()
        qdot_nuc  = empty2D.copy()
        revenue   = empty.copy()
        ppa       = empty.copy()
        annual_e  = empty.copy()
        cap_fac   = empty.copy()
        nuc_cost  = empty.copy()
        fin_cost  = empty.copy()
        TES_CH    = empty.copy()
        TES_DC    = empty.copy()
        iter_log  = empty.copy()
        fail_log  = empty.copy() # 0 = No failure, 1 = SSC caused failure, 
                                 # 2 = Pyomo caused failure, 3 = Pyomo crash occurred
                                 # 4 = SSC low mass flow error, 5 = 
        pyomo_bad_log     = empty.copy()
        pyomo_bad_idx_log = empty.copy()
        pysamdict = {}
        # =============================================================================
        # Operating Modes analysis
        # =============================================================================
        
        from util.DualPlots import DualPlots
        tmp_modes = lambda: None
        DualPlots.set_operating_modes_list(tmp_modes)
        
        op_modes_list = tmp_modes.operating_modes
        
        
        # =============================================================================
        # Double Loop
        # =============================================================================
        
        # starting the time counter
        tic = time.perf_counter()
        
        # TRIPLE LOOP
        for i, th in enumerate(iterator1): #over tshours
            for j, pc in enumerate(iterator2): #over P_cycle rating
                
                # =========================================
                # set up
                # =========================================
                
                # getting indeces
                idx = (i,j)
                mean_idx = (i,j,0)
                stdv_idx = (i,j,1)
                
                # =========================================
                # updating pysam class
                # =========================================
                
                # defining directories
                nuctes = DualPlantTES.DualPlantTES( json_name=json, is_dispatch=dispatch)
        
                # horizons
                nuctes.PySAM_dict['ssc_horizon']   = sscH
                nuctes.ssc_horizon   = sscH * nuctes.u.hr
                nuctes.PySAM_dict['pyomo_horizon'] = pyoH
                nuctes.pyomo_horizon = pyoH * nuctes.u.hr
                
                # saving/updating PYSAM dict to nuctes
                nuctes.dispatch_wrap = nuctes.create_dispatch_wrapper( nuctes.PySAM_dict )
        
                # update SSC dictionary parameters
                nuctes.SSC_dict['P_ref']   = float(pc)
                nuctes.SSC_dict['tshours'] = float(th)
                nuctes.SSC_dict['q_dot_nuclear_des']=q_dot_nuclear_des
                
                """capacity cost for nuclear is 61100000 for 950 MWth nuclear at 450 MWe power cycle
                for solar - assumption is 66/kWyr + 3.5/MWh
                SAM simulations suggest that product is 2288.2 cf 6765 for capacity
                so use capacity cost of 66*(6765+2288.2)/6765 = 88.3
                """
                nuclear_om_fixed = 61100000*q_dot_nuclear_des/950
                solar_om_fixed = 88.3*115000
                if model1:
                    solar_om_fixed=0
                
                
                nuctes.SSC_dict['om_fixed']=[nuclear_om_fixed+solar_om_fixed]
                
                #add in solar paramters
                for param in sol_params:
                    nuctes.SSC_dict[param]=sol[param]
                    
                if case=='zero_drop_cost': #derive trend
                    nuctes.SSC_dict['heliostat_spec_cost']=70
                elif case=='sweep_drop_cost': #to match with Solar PPA to match nuclear PPA
                    nuctes.SSC_dict['heliostat_spec_cost']=43.0
                elif case=='large_drop_cost': 
                    nuctes.SSC_dict['heliostat_spec_cost']=19.1
                    
                #TODO! scale the fixed tower cost if ever want none-115 tower (ref_cost should be OK?)
                
                # update TES costs
                nuctes.SSC_dict["tes_spec_cost"] = tes_spec_cost
                
                #nuctes.SSC_dict['q_dot_rec_des'] = csp_qdot_ref
                
    
                
                nuctes.SSC_dict["plant_spec_cost"]=0.0
                
                # update turbine costs
                
                #Establishes the reference cost of turbine that is bundled into the construction financing cost
                turb_ref_cost=q_dot_nuclear_des/950*400*450
                
                #absolute cost of the turbine
                turb_cost=get_turbine_cost(float(pc))
                turb_premium=turb_cost-turb_ref_cost
                turb_unit_cost=turb_cost/float(pc)
                
                
                #relative cost of the turbine
                nuctes.SSC_dict["bop_spec_cost"]=0.0
                
                #nuclear_spec_cost is believed to be normalized against the turbine, so smear it across the turbine
                nuctes.SSC_dict["nuclear_spec_cost"] = (q_dot_nuclear_des*nuctes.SSC_dict['design_eff']*nuc_spec_cost+turb_cost)/nuctes.SSC_dict['P_ref']
                
                
                # update construction financial costs
                #extra financing of TES should be done against thermal not electric power. This is factor 950/465. NOT IMPLEMENTED for consistency with Model1
                #the impact of this is very small
                base_financing = (nuctes.SSC_dict["construction_financing_cost"]-409000000)+q_dot_nuclear_des/950*409000000
                if model1:
                    base_financing=q_dot_nuclear_des/950*409000000
                extra_financing = fin_yrs * fin_rate * 1000 * (turb_premium + \
                                                        nuctes.SSC_dict['P_ref'] * nuctes.SSC_dict['tshours'] * nuctes.SSC_dict["tes_spec_cost"])
            
                nuctes.SSC_dict["construction_financing_cost"]=base_financing + extra_financing 
        
                # reset design point in Dispatch parameter class
                nuctes.dispatch_wrap.set_design()
                

                # =========================================
                # actually run this simulation
                # =========================================
                
                # print current position in loop
                print(json)
                print("Pref    = {0}".format(nuctes.SSC_dict['P_ref'] ) )
                print("tshours = {0}".format(nuctes.SSC_dict['tshours'] ) )
                
                # =========================================
                # try running simulation, log success
                # =========================================
                
                success = True
                try:
                    # run SSC through PySSC
                    nuctes.run_sim( run_loop=run_loop )
                except Exception as err:
                    print(" Run failed.")
                    success = False
        
                # =========================================
                # log success of simulation
                # =========================================
                
                # extract Plant object
                nt = nuctes.Plant
                
                # log final time index of simulation
                iter_log[idx] = nuctes.slice_ssc_currentH.start if nuctes.slice_ssc_currentH.stop != 8760 \
                                        else nuctes.slice_ssc_currentH.stop
                
                # log success
                success = False if nuctes.slice_ssc_currentH.stop != 8760 else success
                
                # if we ran dispatch optimization, check for failures or crashes 
                if dispatch:
                    
                    # log all of the Pyomo failures (the dictionary stores successes, we negate it with the 'not' command)
                    pyomo_fail = np.array([not nuctes.disp_success[k] for k in nuctes.disp_success.keys()])
                    
                    # number of failures due to Pyomo
                    pyomo_fail_sum = sum(pyomo_fail)
                    
                    # log whether Pyomo failed at all throughout simulation
                    pyomo_bad_log[idx] = bool(pyomo_fail_sum > 0)
                    
                    # log the time index where Pyomo crashed if it did indeed crash
                    if pyomo_fail_sum > 0:
                        pyomo_bad_idx_log[idx] = np.where(pyomo_fail==True)[0][0]
        
                # =========================================
                # Simulation was successful! Log results
                # =========================================
                if success:
                    so = nuctes.SO
                    outputs = DualOutputExtraction(nuctes)
                    
                    # analyzing OP modes
                    op_mode_array       = outputs.op_mode_1
                    op_mode_str_profile = [op_modes_list[ int(x) ] for x in op_mode_array]
                    
                    TES_charging    = np.array(["TES_CH" in x for x in op_mode_str_profile])
                    TES_discharging = np.array(["TES_DC" in x for x in op_mode_str_profile])
                    
                    TES_CH[idx] = int(TES_charging.sum() > 0)
                    TES_DC[idx] = int(TES_discharging.sum() > 0)
                    
                    # log output means
                    gen[mean_idx]       = np.mean( outputs.gen.m )
                    qdot_nuc[mean_idx]  = np.mean( outputs.q_dot_nuc_in.m )
                    defocus[mean_idx]   = np.mean( outputs.defocus )
                    
                    # log output std dev
                    gen[stdv_idx]       = np.std( outputs.gen.m )
                    qdot_nuc[stdv_idx]  = np.std( outputs.q_dot_nuc_in.m )
                    defocus[stdv_idx]   = np.std( outputs.defocus )
                    
                    # log economic outputs
                    revenue[idx]   = np.sum( outputs.gen.to('kW').m * outputs.price )
                    ppa[idx]       = so.Outputs.ppa
                    annual_e[idx]  = (nt.Outputs.annual_energy*u.kWh).to('TWh').m # TWh
                    cap_fac[idx]   = nuctes.capacity_factor.m * 100
                    nuc_cost[idx]  = nt.SystemCosts.nuclear_spec_cost
                    fin_cost[idx]  = so.FinancialParameters.construction_financing_cost
                    
                    print("ppa  = {0} cents/kWh \n\n".format( ppa[idx] ) )
                    print("cost = {0} cents/kWh \n\n".format( fin_cost[idx] ) )
                    print("cap_fac = {0} \n".format(cap_fac[idx]))
                    
                    del outputs, op_mode_array, op_mode_str_profile, so
                
                # =========================================
                # Simulation failed! Log results
                # =========================================
                else:
                    
                    # check if last crash was due to Pyomo
                    if dispatch:
                        # checking the last success log at end of simulation
                        count = int(nuctes.disp_count - 1)
                        disp_success = nuctes.disp_success[count]
                        
                        # Simulation ended because of SSC crash   == 1 
                        # Simulation ended because of Pyomo crash == 2 
                        fail_log[idx] = 1 if disp_success else 2 
                    
                    if hasattr(nuctes,'err_message'):
                        if 'Mass flow rate too low' in nuctes.err_message:
                            fail_log[idx] = 4
                        elif 'stuck' in nuctes.err_message:
                            fail_log[idx] = 5
                        
                    
                    dummy_val = -1
                    # log output means
                    gen[mean_idx]       = dummy_val
                    qdot_nuc[mean_idx]  = dummy_val
                    defocus[mean_idx]   = dummy_val
                    
                    # log output std dev
                    gen[stdv_idx]       = dummy_val
                    qdot_nuc[stdv_idx]  = dummy_val
                    defocus[stdv_idx]   = dummy_val
                    
                    # log economic results
                    revenue[idx]   = dummy_val
                    ppa[idx]       = dummy_val
                    annual_e[idx]  = dummy_val
                    cap_fac[idx]   = dummy_val
                    nuc_cost[idx]  = dummy_val
                    fin_cost[idx]  = dummy_val
                    
                    TES_CH[idx] = dummy_val
                    TES_DC[idx] = dummy_val 
                        
                    # reset the Plant and Grid, prevents memory leak
                del nt
                
                pysamdict = nuctes.PySAM_dict
                
                del nuctes
                    
        # end time counter
        toc = time.perf_counter()        
        
        # print out time elapsed and store it (current record is 5.6 hours)
        time_elapsed = ((toc-tic)*u.s).to('hr')
        print('Made it past the triple loop. Time %.2f hrs' % time_elapsed.m)
        
        
        # =============================================================================
        # Storing Data
        # =============================================================================
        
        # storage dictionary
        Storage = {}
        Storage['time_elapsed']  = time_elapsed
        Storage['tshours']    = tshours
        Storage['p_cycle']    = p_cycle
        Storage['iterator1']  = iterator1
        Storage['iterator2']  = iterator2
        Storage['sscH']       = sscH
        Storage['pyoH']       = pyoH
        Storage['op_modes_list'] = op_modes_list
        Storage['json']       = json
        Storage['dispatch']   = dispatch
        Storage['run_loop']   = run_loop
        Storage['defocus']    = defocus
        Storage['gen']        = gen
        Storage['qdot_nuc']   = qdot_nuc
        Storage['TES_CH']     = TES_CH
        Storage['TES_DC']     = TES_DC
        Storage['iter_log']   = iter_log
        Storage['fail_log']   = fail_log
        Storage['pyomo_bad_log']      = pyomo_bad_log
        Storage['pyomo_bad_idx_log']  = pyomo_bad_idx_log
        Storage['revenue']  = revenue
        Storage['annual_e'] = annual_e
        Storage['ppa']      = ppa
        Storage['cap_fac']  = cap_fac
        Storage['turb_unit_cost']    = turb_unit_cost
        Storage['tes_spec_cost']     = tes_spec_cost
        Storage['fin_yrs']     = fin_yrs
        Storage['fin_rate']    = fin_rate
        Storage['nuclear_spec_cost']            = nuc_cost
        Storage['construction_financing_cost'] = fin_cost
        
        
        # locating output directory
        output_dir = FileMethods.output_dir
        filename = 'paramSweep_varTurbineCost_PySAM__{0}__2022_04__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_[{4},{5}]__PC_[{6},{7}]__{8}.nuctes'.format(
                        json, dispatch, sscH, pyoH, tshours.min(), tshours.max(), p_cycle.min(), p_cycle.max(),case)
        NTPath = os.path.join(output_dir, filename)
        # pickling
        with open(NTPath, 'wb') as f:
            pickle.dump(Storage, f)
