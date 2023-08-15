import pyomo.environ as pyomo
# from pyomo import Parameter, Variable, Set
import matplotlib.pyplot as plt

import math 
import random

# ----------------------------
# prep

N = 96

# Some preliminary sizing calculations
Q_dot_LFR = 950  #MW
W_dot_gen = 475  #MW
eta_cycle_ref = 0.48
cp_salt = 1500 #J/kg-K
DeltaT_ctes = 565-330  #K
DeltaH_csteam = 2125  #kJ/kg  (h(150bar,555C) - h(160bar,300C))
DeltaH_dsteam = 2494  #kJ/kg  (h(4.8bar,x=.98) - h(4.5bar,50C))
DeltaH_dtes = 420 #kJ/kg  (h(5bar,140C) - h(5bar,40C))

m_dot_salt = Q_dot_LFR*1e6 / (cp_salt * DeltaT_ctes)  #kg/s
m_dot_csteam = W_dot_gen/eta_cycle_ref*1e6 / (DeltaH_csteam*1e3)  #kg/s

f_lpt = 0.53  #nominal fraction of steam flow that makes it to LPT stage 3

m_dot_esteam = m_dot_csteam*f_lpt  #nominal extraction flow rate
m_dot_dtes = m_dot_esteam*DeltaH_dsteam/DeltaH_dtes  #nominal flow rate of water into the desal storage tank


ctshours = 5
dtshours = 5

m_ctes_ref = ctshours*m_dot_salt*3600
m_dtes_ref = dtshours*m_dot_dtes*3600


price_schedule = [math.sin(i/24*math.pi ) + 0.5 for i in range(N)]  #electricity price
distillate_schedule = [(math.sin(i/24*math.pi + math.pi/4) + 0.5)*0.2 for i in range(N)]  #distillate price
inflow_schedule = [max(math.sin( i/(N*math.pi) + math.pi/4 )*m_dot_salt/N*random.gauss(1, .25) + m_dot_salt/(N*4), 0) for i in range(N)]  #LFR salt into storage
efficiency_schedule = [0.9-math.cos(i/24*math.pi)*0.1 for i in range(N)]  #ambient efficiency correction

# -------------------------------

model = pyomo.ConcreteModel()

#  ==================================================
#       parameters
#  ==================================================
# model.nt = pyomo.Param(initialize = len(price_schedule), domain = pyomo.Integers)  # number of time steps
model.T = pyomo.Set(initialize = range(1,N+1))  # set of time steps
model.Delta_t = pyomo.Param(initialize = 1)    #[hr]  Time step duration
model.M_dot_LFR = pyomo.Param(model.T, initialize = dict(zip(model.T, inflow_schedule)))    #[kg/s]  Mass flow produced by the LFR at time t
model.eta = pyomo.Param(model.T, initialize = dict(zip(model.T, efficiency_schedule)))   #[]  Ambient temperature efficiency modifier in the power cycle at time t
model.P_elec = pyomo.Param(model.T, initialize = dict(zip(model.T, price_schedule)))    #[$/MWh]  Price at which electricity is sold at time t
model.P_dist = pyomo.Param(model.T, initialize = dict(zip(model.T, distillate_schedule)))   #[$/kg]  Price at which distillate is sold at time t
model.C_c_ramp = pyomo.Param(initialize=0.1)    #[$/Delta MW]  Cost associated with change in electric power production 
model.C_d_ramp = pyomo.Param(initialize=0.2)    #[$/Delta kg/s]  Cost associated with change in distillate production 
model.K_chx  = pyomo.Param(initialize=m_dot_csteam/m_dot_salt)   #[(kg/s)/(kg/s)]  Conversion constant for salt mass flow to steam mass flow in the cycle heat exchanger
model.K_chpt = pyomo.Param(initialize= 0.7*1.1*W_dot_gen)    #[MW/(kg/s)]  Conversion constant for HPT steam mass flow to power in the cycle
model.K_clpt = pyomo.Param(initialize= 0.3*1.1*W_dot_gen)   #[MW/(kg/s)]  Conversion constant for LPT extraction steam mass flow to lost power in the cycle
model.K_cint = pyomo.Param(initialize= -0.1*W_dot_gen)    #[MW]  Constant equation intercept for cycle mass flow to power conversion
model.F_lpt = pyomo.Param(initialize= 0.53)  #[-] Max fraction of steam mass flow available for extraction
model.W_dot_ramp_max = pyomo.Param(initialize=0.1*W_dot_gen)    #[Delta MW/hr]  Maximum allowable power system ramp rate 
# model.V_dot_ramp_max = pyomo.Param(initialize=)    #[Delta m_dot/hr]  Maximum allowable desal system ramp rate 
model.K_dhx = pyomo.Param(initialize=m_dot_esteam/m_dot_dtes)    #[(kg/s)/(kg/s)]  Conversion constant for HPT extraction steam mass flow to desal storage charge mass flow
#model.M_cm_max = pyomo.Param(initialize=m_ctes_ref)    #[kg]  Cycle thermal storage maximum inventory
model.M_cm_min = pyomo.Param(initialize=0)    #[kg]  Cycle thermal storage minimum inventory
#model.M_dm_max = pyomo.Param(initialize=m_dtes_ref)    #[kg]  Desal thermal storage maximum inventory
model.M_dm_min = pyomo.Param(initialize=0)    #[kg]  Desal thermal storage minimum inventory
model.K_d = pyomo.Param(initialize=0.55)    #[(kg/s)/(kg/s)]  Rate of distillate production per unit mass flow into the desal system
model.M_dot_ch_init = pyomo.Param(initialize=0.1*m_ctes_ref)  #[kg] Initial cycle thermal storage inventory
model.M_dot_dh_init  = pyomo.Param(initialize=0.1*m_dtes_ref)  #[kg] Initial desal thermal storage inventory
model.W_dot_max = pyomo.Param(initialize=W_dot_gen)
model.W_dot_min = pyomo.Param(initialize=W_dot_gen*0.25)
model.V_dot_max = pyomo.Param(initialize=m_dot_dtes*model.K_d())
model.V_dot_min = pyomo.Param(initialize=m_dot_dtes*model.K_d()*0)
model.C_s = pyomo.Param(initialize=10) #[$/kg] Cost associated with storage systems

#  ==================================================
#       Variables
#  ==================================================
model.m_ch = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Inventory in cycle-side hot storage at time step $t$
model.m_dot_cs = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow drawn from cycle-side hot storage at time step $t$
model.m_dot_hpt = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow produced in the cycle at the high pressure turbine inlet at time step $t$
model.m_dot_e = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow produced extracted and sent to the desal heat exchanger at time step $t$
model.m_dh = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Inventory in the desal-side hot storage at time step $t$
model.m_dot_ds = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Mass flow consumed by the desal system at time step $t$
model.w_dot = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Electric power produced by the turbine system at time step $t$
model.v_dot = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)   #[]  Distillate produced by the desal system at time step $t$
model.w_dot_Delta_up = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)   #[]  Positive change in electric power produced at time step $t$ relative to the previous time step, $t\geq 2$
model.w_dot_Delta_dn = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)   #[]  Negative change in electric power produced at time step $t$ relative to the previous time step, $t\geq 2$
# model.v_dot_Delta_up = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)   #[]  Positive change in distillate produced at time step $t$ relative to the previous time step, $t\geq 2$
# model.v_dot_Delta_dn = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)   #[]  Negative change in distillate produced at time step $t$ relative to the previous time step, $t\geq 2$
model.M_cm_max = pyomo.Var(model.T, domain=pyomo.NonNegativeReals) #[] Cycle thermal storage maximum inventory
model.M_dm_max = pyomo.Var(model.T, domain=pyomo.NonNegativeReals) #[] Desal thermal storage maximum inventory
#  ==================================================
#   Objective
def objective(model):
    # return sum(model.w_dot[t] for t in model.T)
    return sum([model.P_elec[t]*model.w_dot[t] + model.P_dist[t]*model.v_dot[t] for t in model.T]) \
         - sum([model.C_c_ramp*(model.w_dot_Delta_up[t]+model.w_dot_Delta_dn[t]) for t in (model.T-[1])]) \
         - sum(model.C_s * (model.M_cm_max[t] + model.M_dm_max[t]) for t in model.T) #\
        #  - sum([model.C_d_ramp*(model.v_dot_Delta_up[t]+model.v_dot_Delta_dn[t]) for t in (model.T-[1])])
model.objective = pyomo.Objective(rule=objective, sense=pyomo.maximize)

#  ==================================================
#   Constraints
#  ==================================================

# ------- Power cycle subsystem -------------
def constr_csmass(model,t):
    if t>1: 
        return model.m_ch[t] == (model.M_dot_LFR[t] - model.m_dot_cs[t])*3600*model.Delta_t + model.m_ch[t-1] 
    else:
        return model.m_ch[t] == (model.M_dot_LFR[t] - model.m_dot_cs[t])*3600*model.Delta_t + model.M_dot_ch_init
model.constr_csmass = pyomo.Constraint(model.T, rule=constr_csmass)

# Inventory limits on the cycle storage 
def constr_csmmax(model,t):
    return model.m_ch[t] <= model.M_cm_max[t]
model.constr_csmmax = pyomo.Constraint(model.T, rule=constr_csmmax)
def constr_csmmin(model,t):
    return model.m_ch[t] >= model.M_cm_min
model.constr_csmmin = pyomo.Constraint(model.T, rule=constr_csmmin)

# Mass flow conversion from cycle storage to HPT inlet. \textcolor{red}{Does this need a separate variable, or can we just scale the salt flow?}
def constr_chx(model,t):
    return model.m_dot_hpt[t] == model.K_chx * model.m_dot_cs[t]
model.constr_chx = pyomo.Constraint(model.T, rule=constr_chx)

# 

# Conversion of mass flow to power in the cycle
def constr_cmass(model, t):
    return model.w_dot[t] == model.eta[t] * (model.K_chpt * model.m_dot_hpt[t] - model.K_clpt * model.m_dot_e[t] + model.K_cint)
model.constr_cmass = pyomo.Constraint(model.T, rule=constr_cmass)

# Mass flow limit on the extraction
def constr_mext(model, t):
    return model.m_dot_e[t] <= model.F_lpt * model.m_dot_hpt[t]
model.constr_mext = pyomo.Constraint(model.T, rule=constr_mext)

# Range limits on power
def constr_wbounds_up(model, t):
    return model.w_dot[t] <= model.W_dot_max
model.constr_wbounds_up = pyomo.Constraint(model.T, rule=constr_wbounds_up)
def constr_wbounds_dn(model, t):
    return model.w_dot[t] >= model.W_dot_min
model.constr_wbounds_dn = pyomo.Constraint(model.T, rule=constr_wbounds_dn)

# Power ramping tracking
def constr_wtrack_up(model,t):
    return model.w_dot_Delta_up[t] >= model.w_dot[t] - model.w_dot[t-1]
model.constr_wtrack_up = pyomo.Constraint(model.T-[1], rule=constr_wtrack_up)
def constr_wtrack_dn(model,t):
    return model.w_dot_Delta_dn[t] >= model.w_dot[t-1] - model.w_dot[t]
model.constr_wtrack_dn = pyomo.Constraint(model.T-[1], rule=constr_wtrack_dn)

# Power ramping limits
def constr_cramp_up(model,t):
    return model.w_dot_Delta_up[t] <= model.W_dot_ramp_max
model.constr_cramp_up = pyomo.Constraint(model.T-[1], rule=constr_cramp_up)
def constr_cramp_dn(model,t):
    return model.w_dot_Delta_dn[t] <= model.W_dot_ramp_max
model.constr_cramp_dn = pyomo.Constraint(model.T-[1], rule=constr_cramp_dn)



# ------- Desal subsystem -------------
# Mass balance on the desal storage
def constr_dsmass(model,t):
    if t>1: 
        return model.m_dh[t] == model.K_dhx*model.m_dot_e[t] + model.m_dh[t-1] - model.m_dot_ds[t]
    else:
        return model.m_dh[t] == model.K_dhx*model.m_dot_e[t] + model.M_dot_dh_init - model.m_dot_ds[t]
model.constr_dsmass = pyomo.Constraint(model.T, rule=constr_dsmass)

# Inventory limits on the desal storage 
def constr_dsmmax(model,t):
    return model.m_dh[t] <= model.M_dm_max[t]
model.constr_dsmmax = pyomo.Constraint(model.T, rule=constr_dsmmax)
def constr_dsmmin(model,t):
    return model.m_dh[t] >= model.M_dm_min
model.constr_dsmmin = pyomo.Constraint(model.T, rule=constr_dsmmin)


# Mass flow conversion from desal storage to distillate production 
def constr_distillate(model,t):
    return model.v_dot[t] == model.K_d * model.m_dot_ds[t]
model.constr_distillate = pyomo.Constraint(model.T, rule=constr_distillate)

# Range limits on desal system
def constr_vbounds_up(model, t):
    return model.v_dot[t] <= model.V_dot_max
model.constr_vbounds_up = pyomo.Constraint(model.T, rule=constr_vbounds_up)
def constr_vbounds_dn(model, t):
    return model.v_dot[t] >= model.V_dot_min
model.constr_vbounds_dn = pyomo.Constraint(model.T, rule=constr_vbounds_dn)

# Desal ramping tracking
# \dot{v}^{\Delta +}_t \geq \dot{v}_t - \dot{v}_{t-1}\ \ &\forall t \in \mathcal{T}\ :\ t \geq 2 \\ 
# \dot{v}^{\Delta -}_t \geq \dot{v}_{t-1} - \dot{v}_{t}\ \ &\forall t \in \mathcal{T}\ :\ t \geq 2 \\ 
# \dot{v}^{\Delta -}_t \geq 0,\ \dot{v}^{\Delta +}_t \geq 0 \ \ &\forall t \in \mathcal{T}  

# Desal ramping limits
# \dot{v}^{\Delta +} \leq \dot{V}^{ramp,max}\ \ &\forall t \in \mathcal{T} \\
# \dot{v}^{\Delta -} \leq \dot{V}^{ramp,max}\ \ &\forall t \in \mathcal{T} 
# \label{eq:vramp}

# ------------------


solver = pyomo.SolverFactory('gurobi')
results = solver.solve(model, keepfiles = True, logfile = 'solve.log')
# print(model.display())
print(results)
outlabs = ['time', 'P_elec','P_dist','m_ch', 'm_dh', 'w_dot', 'v_dot', 'm_dot_cs','m_dot_hpt','m_dot_e','m_dot_ds', 'M_cm_max', 'M_dm_max']
print( ('{:6s}\t'+ '\t'.join(['{:10s}' for i in range(len(outlabs)-1)])).format(*outlabs) )
for t in model.T:
    outs = [
        t, 
        model.P_elec[t], 
        model.P_dist[t], 
        model.m_ch[t]()*1e-3, 
        model.m_dh[t]()*1e-3, 
        model.w_dot[t](), 
        model.v_dot[t](),
        model.m_dot_cs[t](),
        model.m_dot_hpt[t](),
        model.m_dot_e[t](),
        model.m_dot_ds[t](),
        model.M_cm_max[t]()*1e-3
        model.M_dm_max[t]()*1e-3
        ]
    fstr = '{:d}\t' + '\t'.join(['{:f}' for i in range(len(outs)-1)])
    print(fstr.format(*outs))

