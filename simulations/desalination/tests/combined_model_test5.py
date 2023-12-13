#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 21:26:00 2023

@author: elizabethkeith
"""

"For this test case, the price at which water is sold at is set to 0, as is the penalty associated \
    with ramping the NPP, maximum output from NPP is twice the nominal output"

import pyomo.environ as pyomo
# from pyomo import Parameter, Variable, Set
import pandas as pd
import numpy as np

# ----------------------------
# prep

N                         = 8760

# Some preliminary sizing calculations
Q_dot_LFR                 = 950  #MW
W_dot_gen                 = 494 #MW
eta_cycle_ref             = 0.52
cp_salt                   = 1.5 #kJ/kg-K
DeltaT_ctes               = 565-330  #K
DeltaH_csteam             = 2126  #kJ/kg  (h(150bar,555C) - h(150bar,300C))
DeltaH_dsteam             = 2494  #kJ/kg  (h(4.8bar,x=.98) - h(4.5bar,50C))  # steam extracted from turbine
DeltaH_dtes               = 421.3 #kJ/kg  (h(5bar,140C) - h(5bar,40C))       # water through HX to desal hot storage 

# calculations for molten salt TES
cost_molten_salt_TES_kWh  = 10   # [$/kWh-th]
#cost_molten_salt_TES_kg   = (cost_molten_salt_TES_kWh/3600)*cp_salt*DeltaT_ctes  # [$/kg]

# calculations for desalination TES
cost_desal_TES_kWh        = 1  # [$/KWh-th], scaled from German TES
#cost_desal_TES_kg         = (cost_desal_TES_kWh/3600)*DeltaH_dtes # [$/kg]

# esalinated water prices
desal_price_acre_foot     = 2637.10 # [$/acre-foot]  price of desalinated water in 2021 in San Diego County from:(https://www.sdcwa.org/wp-content/uploads/2020/11/desal-carlsbad-fs.pdf)
cubic_meter_per_acre_foot = 1233.48 # [m^3/acre-foot]
specific_volume_H2O       = 0.001 # [m^3/kg]
desal_price_per_kg        = 1


m_dot_salt                = Q_dot_LFR*1e3 / (cp_salt * DeltaT_ctes)  #kg/s
m_dot_csteam              = W_dot_gen/eta_cycle_ref*1e6 / (DeltaH_csteam*1e3)  #kg/s
#m_dot_csteam_ovsz         = 2*W_dot_gen/eta_cycle_ref*1e6 / (DeltaH_csteam*1e3)

f_lpt                     = 0.53  #nominal fraction of steam flow that makes it to LPT stage 3

m_dot_esteam              = m_dot_csteam*f_lpt  #nominal extraction flow rate
m_dot_dtes                = m_dot_esteam*DeltaH_dsteam/DeltaH_dtes  #nominal flow rate of water into the desal storage tank


distillate_price_schedule = [desal_price_per_kg for i in range(N)]  #distillate price in San Diego ($/kg) (https://www.sdcwa.org/wp-content/uploads/2020/11/desal-carlsbad-fs.pdf)
inflow_schedule           = [m_dot_salt for i in range(N)] #LFR salt into storage
temperature_schedule      = pd.read_csv("./Diablo_Canyon_Temps_2021.csv", usecols=["Temperature"]) #temperatures in Diablo Canyon for TMY
temperature_schedule      = temperature_schedule["Temperature"].values

efficiency_schedule       = []    #efficiency schedule as function of ambient temperatures
for temp in temperature_schedule:
    efficiency = (-0.0036*temp + 1.2161) - 0.073
    efficiency_schedule.append(efficiency)

electric_price_schedule   = pd.read_csv("../../electric_price_schedule.csv")
electric_price_schedule   = electric_price_schedule["LMP"].values

# # -------------------------------

model = pyomo.ConcreteModel()

# #  ==================================================
# #       parameters
# #  ==================================================
model.nt             = pyomo.Param(initialize = len(inflow_schedule), domain = pyomo.Integers)  # number of time steps
model.T              = pyomo.Set(initialize = range(1,N+1))  # set of time steps
model.Delta_t        = pyomo.Param(initialize = 1)    #[hr]  Time step duration
model.M_dot_LFR      = pyomo.Param(model.T, initialize = dict(zip(model.T, inflow_schedule)))    #[kg/s]  Mass flow produced by the LFR at time t
model.eta            = pyomo.Param(model.T, initialize = dict(zip(model.T, efficiency_schedule)))   #[]  Ambient temperature efficiency modifier in the power cycle at time t
model.P_elec         = pyomo.Param(model.T, initialize = dict(zip(model.T, electric_price_schedule)))    #[$/MWh]  Price at which electricity is sold at time t
model.P_dist         = pyomo.Param(model.T, initialize = dict(zip(model.T, distillate_price_schedule)))   #[$/kg]  Price at which distillate is sold at time t
model.C_c_ramp       = pyomo.Param(initialize=0)    #[$/Delta MW]  Cost associated with change in electric power production, from Gabriel's paper
model.C_d_ramp       = pyomo.Param(initialize=0.2)    #[$/Delta kg/s]  Cost associated with change in distillate production 
model.K_chx          = pyomo.Param(initialize=m_dot_csteam/m_dot_salt)   #[(kg/s)/(kg/s)]  Conversion constant for salt mass flow to steam mass flow in the cycle heat exchanger
#model.K_chpt         = pyomo.Param(initialize=W_dot_gen/m_dot_csteam)    #[MW/(kg/s)]  Conversion constant for HPT steam mass flow to power in the cycle
#model.K_clpt         = pyomo.Param(initialize=W_dot_gen/m_dot_esteam)   #[MW/(kg/s)]  Conversion constant for LPT extraction steam mass flow to lost power in the cycle
model.K_cint         = pyomo.Param(initialize=-0.073*W_dot_gen)    #[MW]  Constant equation intercept for cycle mass flow to power conversion
model.F_lpt          = pyomo.Param(initialize=0.53)  #[-] Max fraction of steam mass flow available for extraction
model.W_dot_ramp_max = pyomo.Param(initialize=0.1*W_dot_gen)    #[Delta MW/hr]  Maximum allowable power system ramp rate 
model.K_dhx          = pyomo.Param(initialize=m_dot_dtes/m_dot_esteam)    #[(kg/s)/(kg/s)]  Conversion constant for HPT extraction steam mass flow to desal storage charge mass flow
model.M_cm_min       = pyomo.Param(initialize=0)    #[kg]  Cycle thermal storage minimum inventory
model.M_dm_min       = pyomo.Param(initialize=0)    #[kg]  Desal thermal storage minimum inventory
model.K_d            = pyomo.Param(initialize=0.55)    #[(kg/s)/(kg/s)]  Rate of distillate production per unit mass flow into the desal system
model.M_dot_ch_init  = pyomo.Param(initialize=0.0)  #[kg] Initial cycle thermal storage inventory
model.M_dot_dh_init  = pyomo.Param(initialize=0.0)  #[kg] Initial desal thermal storage inventory
model.W_dot_max      = pyomo.Param(initialize=W_dot_gen*2) #[MW]   
model.W_dot_min      = pyomo.Param(initialize=W_dot_gen*0.25)    #[MW]
# model.V_dot_max      = pyomo.Param(initialize=m_dot_dtes*model.K_d()) #[kg/s]
model.V_dot_min      = pyomo.Param(initialize=m_dot_dtes*model.K_d()*0) #[kg/s]
model.C_cs           = pyomo.Param(initialize=cost_molten_salt_TES_kWh) #[$/kWh] Cost associated with cycle thermal energy storage
model.C_ds           = pyomo.Param(initialize=cost_desal_TES_kWh) #[$/kWh] Cost associated with desal thermal energy storage
#model.m_dot_csteam  = pyomo.Param(initialize=m_dot_csteam) #[kg/s] Steam flow to turbine
#model.m_dot_cs      = pyomo.Param(initialize=m_dot_salt)
#model.m_dot_hpt      = pyomo.Param(initialize=m_dot_csteam)
#model.M_cm_max       = pyomo.Param(initialize=2E6) #2000 MWh of cycle-side storage
model.C_desal        = pyomo.Param(initialize=30E6)
#  ==================================================
#       Variables
#  ==================================================
model.m_ch           = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Inventory in cycle-side hot storage at time step $t$
model.m_dot_cs       = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow drawn from cycle-side hot storage at time step $t$
model.m_dot_hpt      = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow produced in the cycle at the high pressure turbine inlet at time step $t$
model.m_dot_e        = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow produced extracted and sent to the desal heat exchanger at time step $t$
model.m_dh           = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Inventory in the desal-side hot storage at time step $t$
model.m_dot_ds       = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow consumed by the desal system at time step $t$
model.w_dot          = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Electric power produced by the turbine system at time step $t$
model.v_dot          = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Distillate produced by the desal system at time step $t$
model.w_dot_Delta_up = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)   #[]  Positive change in electric power produced at time step $t$ relative to the previous time step, $t\geq 2$
model.w_dot_Delta_dn = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)   #[]  Negative change in electric power produced at time step $t$ relative to the previous time step, $t\geq 2$
model.M_cm_max       = pyomo.Var(domain=pyomo.NonNegativeReals) #[kWh] Cycle thermal storage maximum inventory
model.M_dm_max       = pyomo.Var(domain=pyomo.NonNegativeReals) #[kWh] Desal thermal storage maximum inventory
model.V_dot_max      = pyomo.Var(domain=pyomo.NonNegativeReals) #[] Size of desalination system
#model.C_desal        = pyomo.Var(domain=pyomo.NonNegativeReals) #[] Cost of desalination system
#  ==================================================
#   Objective
def objective(model):
    # return sum(model.w_dot[t] for t in model.T)
    return sum([model.P_elec[t]*model.w_dot[t]*model.Delta_t + model.P_dist[t]*model.v_dot[t]*3600*model.Delta_t for t in model.T]) \
        - sum([model.C_c_ramp*(model.w_dot_Delta_up[t]+model.w_dot_Delta_dn[t]) for t in (model.T-[1])]) \
        - (model.C_cs*model.M_cm_max)/30 - (model.C_ds*model.M_dm_max)/30 - model.C_desal/30
model.objective = pyomo.Objective(rule=objective, sense=pyomo.maximize)

#  ==================================================
#   Constraints
#  ==================================================

def constr_Vdotmax(model,t):
    return model.V_dot_max <= model.v_dot[t] 
model.constr_Vdotmax = pyomo.Constraint(model.T, rule=constr_Vdotmax)

def constr_C_desal1(model):
    return model.C_desal >= -0.003*model.V_dot_max + 1300.1
model.constr_C_desal1 = pyomo.Constraint(rule=constr_C_desal1)

def constr_C_desal2(model):
    return model.C_desal >= -0.0218*model.V_dot_max + 1648.8
model.constr_C_desal2 = pyomo.Constraint(rule=constr_C_desal2)



# ------- Power cycle subsystem ------------- 
def constr_csmass(model,t):
    if t>1: 
        return model.m_ch[t] == (model.M_dot_LFR[t] - model.m_dot_cs[t]) * 3600 * model.Delta_t + model.m_ch[t-1]
    else:
        return model.m_ch[t] == (model.M_dot_LFR[t] - model.m_dot_cs[t]) * 3600 * model.Delta_t + model.M_dot_ch_init
model.constr_csmass = pyomo.Constraint(model.T, rule=constr_csmass)

#Inventory limits on the cycle storage 
def constr_csmmax(model,t):
    return (model.m_ch[t] * cp_salt * DeltaT_ctes)/3600 <= model.M_cm_max
model.constr_csmmax = pyomo.Constraint(model.T, rule=constr_csmmax)
def constr_csmmin(model,t):
    return (model.m_ch[t] * cp_salt * DeltaT_ctes)/3600 >= model.M_cm_min
model.constr_csmmin = pyomo.Constraint(model.T, rule=constr_csmmin)


#Mass flow conversion from cycle storage to HPT inlet
def constr_chx(model,t):
    return model.m_dot_hpt[t] == (model.K_chx * model.m_dot_cs[t])
model.constr_chx = pyomo.Constraint(model.T, rule=constr_chx)

def constr_hptmax(model,t):
    return model.m_dot_hpt[t] <= 2*W_dot_gen/eta_cycle_ref*1e6 / (DeltaH_csteam*1e3)
model.constr_hptmax = pyomo.Constraint(model.T, rule=constr_hptmax)


# Conversion of mass flow to power in the cycle
def constr_cmass(model, t):
    return model.w_dot[t] == model.eta[t] * W_dot_gen * (0.7 * (model.m_dot_hpt[t]/m_dot_csteam) + 0.3*(0.62*(1 - model.m_dot_e[t]/m_dot_esteam)+0.38))
model.constr_cmass = pyomo.Constraint(model.T, rule=constr_cmass)
                                                         

# Mass flow limit on the extraction
def constr_mext(model, t):
    return model.m_dot_e[t] <= model.F_lpt * model.m_dot_hpt[t]
model.constr_mext = pyomo.Constraint(model.T, rule=constr_mext)

# Range limits on power
def constr_wbounds_up(model, t):
    return model.w_dot[t] <= model.W_dot_max * model.eta[t]
model.constr_wbounds_up = pyomo.Constraint(model.T, rule=constr_wbounds_up)
def constr_wbounds_dn(model, t):
    return model.w_dot[t] >= model.W_dot_min
model.constr_wbounds_dn = pyomo.Constraint(model.T, rule=constr_wbounds_dn)

# Power ramping tracking
def constr_wtrack_up(model, t):
    return model.w_dot_Delta_up[t] >= model.w_dot[t] - model.w_dot[t-1]
model.constr_wtrack_up = pyomo.Constraint(model.T-[1], rule=constr_wtrack_up)
def constr_wtrack_dn(model, t):
    return model.w_dot_Delta_dn[t] >= model.w_dot[t-1] - model.w_dot[t]
model.constr_wtrack_dn = pyomo.Constraint(model.T-[1], rule=constr_wtrack_dn)

# Power ramping limits
def constr_cramp_up(model,t):
    return model.w_dot_Delta_up[t] <= model.W_dot_ramp_max*model.Delta_t
model.constr_cramp_up = pyomo.Constraint(model.T-[1], rule=constr_cramp_up)
def constr_cramp_dn(model,t):
    return model.w_dot_Delta_dn[t] <= model.W_dot_ramp_max*model.Delta_t
model.constr_cramp_dn = pyomo.Constraint(model.T-[1], rule=constr_cramp_dn)



# ------- Desal subsystem -------------
# Mass balance on the desal storage
def constr_dsmass(model,t):
    if t>1: 
        return model.m_dh[t] == (model.K_dhx*model.m_dot_e[t]-model.m_dot_ds[t])*3600*model.Delta_t + model.m_dh[t-1]
    else:
        return model.m_dh[t] == (model.K_dhx*model.m_dot_e[t]-model.m_dot_ds[t])*3600*model.Delta_t + model.M_dot_dh_init
model.constr_dsmass = pyomo.Constraint(model.T, rule=constr_dsmass)

# Inventory limits on the desal storage 
def constr_dsmmax(model,t):
    return model.m_dh[t] <= (model.M_dm_max*3600)/DeltaH_dtes
model.constr_dsmmax = pyomo.Constraint(model.T, rule=constr_dsmmax)
def constr_dsmmin(model,t):
    return model.m_dh[t] >= model.M_dm_min
model.constr_dsmmin = pyomo.Constraint(model.T, rule=constr_dsmmin)


# Mass flow conversion from desal storage to distillate production 
def constr_distillate(model,t):
    return model.v_dot[t] == model.K_d * model.m_dot_ds[t]
model.constr_distillate = pyomo.Constraint(model.T, rule=constr_distillate)

# # Range limits on desal system
# def constr_vbounds_up(model, t):
#     return model.v_dot[t] <= model.V_dot_max
# model.constr_vbounds_up = pyomo.Constraint(model.T, rule=constr_vbounds_up)
def constr_vbounds_dn(model, t):
    return model.v_dot[t] >= model.V_dot_min
model.constr_vbounds_dn = pyomo.Constraint(model.T, rule=constr_vbounds_dn)

# ------------------


solver = pyomo.SolverFactory('gurobi')
results = solver.solve(model, keepfiles = True, logfile = 'solve.log')
# print(model.display())

outlabs = ['t', 'P_elec','eta','m_ch', 'm_dh', 'w_dot', 'v_dot', 'm_dot_cs','m_dot_hpt','m_dot_e','m_dot_ds', 'M_cm_max', 'M_dm_max', 'V_dot_max']
data_out = np.zeros((8760,14))
for t in model.T:
    outs = [
        t, 
        model.P_elec[t],
        model.eta[t],
        model.m_ch[t](), 
        model.m_dh[t](), 
        model.w_dot[t](),
        model.v_dot[t](),
        model.m_dot_cs[t](),
        model.m_dot_hpt[t](),
        model.m_dot_e[t](),
        model.m_dot_ds[t](),
        model.M_cm_max(),
        model.M_dm_max(),
        model.V_dot_max(),

        ]

    data_out[t-1] = np.array(outs)
df_out = pd.DataFrame(data_out, columns=outlabs)
df_out.to_csv('results5.csv', header=True)



