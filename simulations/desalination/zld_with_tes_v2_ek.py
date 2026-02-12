import pyomo.environ as pyo
import pandas as pd
import os


## for all of these cases, the electricity price is set to $0.12/kWhe
# BASE CASE: SW to MED turned on

# TEST 1: Decrease RO electricity requirement
# Result: When decreased to 3.0 kWe/m3, RO did not turn on; however, decreasing to 2.5 kWe/m3 the same results happen as steady state
#         model where RO turns on and goes RO->MED


# TEST 2: 10x MED electricity use
# Results: Does not initially turn on only RO like in base case; but when I change the multiplier on water to 1.2 and greater it does exhibit the same behavior


# TEST 3: Water is worth nothing
# Result: No active links, same as steady-state model


# TEST 4: Lithium price is x100
# Result: SW->RO->MED->ED (different from steady-state where it is SW->MED->ED)


# TEST 5: Lithium price is x100 and MED electricity is x10
# Result: SW->RO->ED (different from steady-state where it is SW->RO->MED->ED)


# TEST 6: Lithium price is x100 and MED electricity is x100
# Result: SW->RO->ED (same as steady-state model)


# TEST 7: Salt price 100$/kg
# Result: SW->RO->MED->CRY (same as steady-state model)


# TEST 8: Salt price $100/kg and lithium price 100x base cost
# Result: SW->RO->MED->ED->CRY (same as base model)


# TEST 9: 10x MED electricity, RO at 4.6 kWe/m3
# Result: No desalination turned on (same as steady-state)


# TEST 10: Oversized turbine (x2), RO elec x100, cry elec x100 (so those aren't built)
# Result: Sawtooth harging and discharging of tes systems, constant med output
    
    

from zld_params import (
    add_sets,
    add_time_set,
    add_time_series_params,
    add_thermo_params,
    add_process_params,
    add_economic_params,
    add_zld_params
)

Schedule_effic = pd.read_csv("sched_eff.csv", usecols=["Efficiency"]).to_numpy().flatten()
Schedule_elec  = pd.read_csv("DIABLOCN_2_N001_PRICES.csv", usecols=["LMP"]).to_numpy().flatten()


model = pyo.ConcreteModel()

N = 168

add_sets(model)
add_time_set(model, N)
add_time_series_params(model, sched_price_elec=Schedule_elec, sched_cyc_eff=Schedule_effic)
add_thermo_params(model, oversize_factor=1.0)
add_process_params(model, process_scale={'RO_elec':1.0, 'CRY_elec':1})
add_economic_params(model, price_scale={'water'  : 1.0,
                                        'lithium': 1.0}, cost_scale={})
add_zld_params(model, yield_scale=1.0)



# variables 

model.m_hot_tes               = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kg,            inventory in high-temp hot storage 
model.m_dot_hot_tes           = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kg/s,          mass flow out of high-temp hot storage
model.m_dot_hpt               = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kg/s,          mass flow at the high pressure turbine inlet 
model.m_dot_extract           = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kg/s,          mass flow of steam extracted from power cycle
model.m_cold_tes              = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kg,            inventory in low-temp hot storage
model.m_dot_cold_tes          = pyo.Var(model.T, model.processes, domain=pyo.NonNegativeReals)              # units: kg/s,          mass flow out of low-temperature storage to med
model.w_dot_gen               = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kW_e,          electricity produced by power cycle
model.m_hot_tes_max           = pyo.Var(domain=pyo.NonNegativeReals)                                        # units: kWh_th,        high-temp thermal storage maximum inventory
model.m_cold_tes_max          = pyo.Var(domain=pyo.NonNegativeReals)                                        # units: kWh_th,        low-temp thermal storage maximum inventory
model.m_dot_hpt_Delta_up      = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kg/s,          upward change in hpt mass flow from one time step to the next 
model.m_dot_hpt_Delta_dn      = pyo.Var(model.T, domain=pyo.NonNegativeReals)                               # units: kg/s,          downward change in hpt mass flow from one time step to the next
model.y_link_active           = pyo.Var(model.links, within=pyo.Binary)                                     # units: --,            1 if there is flow from process q to process p; 0 otherwise
model.y_process_active        = pyo.Var(model.processes, within=pyo.Binary)                                 # units: --,            1 if process p is active; 0 otherwise
model.v_dot_in                = pyo.Var(model.T, model.processes, within=pyo.NonNegativeReals)              # units: m^3/hr,        volumetric flow rate entering each process
model.v_dot_conc              = pyo.Var(model.T, model.processes, within=pyo.NonNegativeReals)              # units: m^3/hr,        volumetric flow rate of concentrated stream leaving each process
model.v_dot_dil               = pyo.Var(model.T, model.processes, within=pyo.NonNegativeReals)              # units: m^3/hr,        volumetric flow rate of diluted stream leaving each process
model.v_dot_link              = pyo.Var(model.T, model.links, within=pyo.NonNegativeReals)                  # units: m^3/hr,        volumetric flow rate from an upstream process to a downstream process
model.elec_used               = pyo.Var(model.T, model.processes, within=pyo.NonNegativeReals)              # units: kW_e,          electricity used by each process
model.elec_sold               = pyo.Var(model.T, within=pyo.NonNegativeReals)                               # units: kW_e,          electricity sold to the grid from each power cycle configuration   
model.q_generated             = pyo.Var(model.T, within=pyo.NonNegativeReals)                               # units: kWh_th,        heat generated by each power cycle from steam extractions
model.q_used                  = pyo.Var(model.T, model.processes, within=pyo.NonNegativeReals)              # units: kWh_th,        heat used by each process
model.f_ion_in                = pyo.Var(model.T, model.processes, model.ions, within=pyo.NonNegativeReals)  # units: kg/h,          flow of each ion into each process
model.f_ion_conc              = pyo.Var(model.T, model.processes, model.ions, within=pyo.NonNegativeReals)  # units: kg/h,          ion flow in concentrated stream leaving process p
model.f_ion_dil               = pyo.Var(model.T, model.processes, model.ions, within=pyo.NonNegativeReals)  # units: kg/h,          ion flow in diluted stream leaving process p
model.f_ion_link              = pyo.Var(model.T, model.links, model.ions, within=pyo.NonNegativeReals)      # units: kg/h,          ion flow between two processes
model.li_recovered_ed         = pyo.Var(model.T, within=pyo.NonNegativeReals)                               # units: kg/h,          total Li recovered
model.v_dot_cap               = pyo.Var(model.processes, within=pyo.NonNegativeReals)                       # units: m^3/day,       capacity of each process





# convert steam to power in the power cycle
def steam_to_power(model, t):
    return model.w_dot_gen[t] == (
        model.Eff_cycle[t]*model.W_dot_ref*
        (model.Slope_power_curve*model.m_dot_hpt[t]/model.M_dot_hx1_steam_nom
            - model.F_lpt_post*model.m_dot_extract[t]/model.M_dot_hx1_steam_nom
            - model.Interc_power_curve))
model.steam_to_power = pyo.Constraint(model.T, rule=steam_to_power)



# limits extraction from power cycle
def steam_extraction_limit(model, t):
    return model.m_dot_extract[t] <= model.F_ext * model.m_dot_hpt[t]
model.steam_extraction_limit = pyo.Constraint(model.T, rule=steam_extraction_limit)



# limits flow to hpt depending on turbine size
def bound_hpt_flow(model, t):
    return model.m_dot_hpt[t] <= model.M_dot_hpt_max
model.bound_hpt_flow = pyo.Constraint(model.T, rule=bound_hpt_flow)



# relates salt mass flow out of hot storage to steam flow in power cycle
def hx1_conversion(model,t):
    return model.m_dot_hpt[t] == model.K_hx1 * model.m_dot_hot_tes[t]
model.hx1_conversion = pyo.Constraint(model.T, rule=hx1_conversion)



# limits hot tes inventory
def hot_tes_limit(model,t):
    return model.m_hot_tes[t] <= model.m_hot_tes_max
model.hot_tes_limit = pyo.Constraint(model.T, rule=hot_tes_limit)



# limits cold tes inventory
def cold_tes_limit(m, t):
    return m.m_cold_tes[t] <= m.m_cold_tes_max
model.cold_tes_limit = pyo.Constraint(model.T, rule=cold_tes_limit)



def prev_circular(T, t):
    return T.prev(t) if t != T.first() else T.last()

# tracks power ramping upwards 
def hpt_track_up(model, t):
    t_prev = prev_circular(model.T,t)
    return model.m_dot_hpt_Delta_up[t] >= model.m_dot_hpt[t] - model.m_dot_hpt[t_prev]
model.hpt_track_up = pyo.Constraint(model.T, rule=hpt_track_up)


# tracks power ramping downwards
def hpt_track_down(model, t):
    t_prev = prev_circular(model.T,t)
    return model.m_dot_hpt_Delta_dn[t] >= model.m_dot_hpt[t_prev] - model.m_dot_hpt[t]
model.hpt_track_down = pyo.Constraint(model.T, rule=hpt_track_down)



# limits power ramping upwards based on power ramp limit
def hpt_ramp_up(model,t):
    return model.m_dot_hpt_Delta_up[t] <= model.W_dot_ramp_max/model.K_power*model.Delta_t
model.hpt_ramp_up = pyo.Constraint(model.T, rule=hpt_ramp_up)



# limits power ramping downwards based on power ramp limit
def hpt_ramp_down(model,t):
    return model.m_dot_hpt_Delta_dn[t] <= model.W_dot_ramp_max/model.K_power*model.Delta_t
model.hpt_ramp_down = pyo.Constraint(model.T, rule=hpt_ramp_down)  



# conserves mass in hot tes
def hot_tes_inventory(model, t):  
    t_prev = prev_circular(model.T,t)
    return model.m_hot_tes[t] == (model.M_dot_salt[t] - model.m_dot_hot_tes[t]) * model.Convert_time * model.Delta_t + model.m_hot_tes[t_prev]
model.hot_tes_inventory = pyo.Constraint(model.T, rule=hot_tes_inventory)



# conserves mass in cold tes low-temp tes inventory tracking between time steps
def cold_tes_inventory(model, t):
    t_prev = prev_circular(model.T,t)
    return model.m_cold_tes[t] == model.m_cold_tes[t_prev] + (
            model.K_hx2 * model.m_dot_extract[t] 
            - sum(model.m_dot_cold_tes[t, p] for p in ['MED', 'CRY'])) * model.Delta_t
model.cold_tes_inventory = pyo.Constraint(model.T, rule=cold_tes_inventory)



# relates electricity used by each process to per unit requirements
def electric_used_by_process(model, t, p):
    return model.elec_used[t,p] == model.Elec_required[p] * model.v_dot_in[t,p]
model.electric_used_by_process = pyo.Constraint(model.T, model.processes, rule=electric_used_by_process)



# balances electricity 
def electric_balance(model, t):
    return model.elec_sold[t] == model.w_dot_gen[t] - sum(model.elec_used[t, p] for p in model.processes)
model.electric_balance = pyo.Constraint(model.T, rule=electric_balance)



# routes cold storage to med or cry
def cold_tes_balance(model, t, p):
    return model.m_dot_cold_tes[t, p] * model.Delta_H_hx2_water * model.Delta_t == model.q_used[t, p]
model.cold_tes_balance = pyo.Constraint(model.T, model.processes, rule=cold_tes_balance)



# relates heat used by each process to per unit requirements
def heat_used_by_process(model, t, p):
    return model.q_used[t, p] == model.Q_required[p] * model.v_dot_in[t, p]
model.heat_used_by_process = pyo.Constraint(model.T, model.processes, rule=heat_used_by_process)



# conserves ions through system
def ion_conservation(m, t, p, i):
    return m.f_ion_in[t,p,i] == m.f_ion_conc[t,p,i] + m.f_ion_dil[t,p,i]
model.ion_conservation = pyo.Constraint(model.T, model.processes, model.ions, rule=ion_conservation)



# relates ions into a process as sum of ions on upstream links
def ion_mass_in(m,t,p,i):
    return m.f_ion_in[t,p,i] == sum(m.f_ion_link[t,(q,p),i]for (q,r) in m.links if r == p)
model.ion_mass_in = pyo.Constraint(model.T, model.processes, model.ions, rule=ion_mass_in)



# relates ion flow on link between seawater and process to seawater concentration
def sw_ion_link(m,t,p,i):
    if ('SW',p) not in m.links:
        return pyo.Constraint.Skip
    return m.f_ion_link[t,('SW',p),i] == m.v_dot_link[t,('SW',p)] * m.Conc_sw[i]
model.sw_ion_link = pyo.Constraint(model.T, model.processes, model.ions, rule=sw_ion_link)



# constrains flow only on active links
def link_gate_all(m, t, p, q):
    return m.v_dot_link[t,(p,q)] <= m.V_dot_max * m.y_link_active[(p,q)]
model.link_gate_all = pyo.Constraint(model.T, model.links, rule=link_gate_all)



# constrains flow to one downstream process
def one_outgoing(m, t, p):
    if p == 'SW':
        return pyo.Constraint.Skip

    outs = [q for (p2, q) in m.links if p2 == p]
    if not outs:
        return pyo.Constraint.Skip

    return sum(m.y_link_active[(p, q)] for q in outs) <= 1

model.one_outgoing = pyo.Constraint(model.T, model.processes, rule=one_outgoing)



# conserves volumetric flow through each process
def process_flow_balance(model, t, p):
    return model.v_dot_in[t, p] == model.v_dot_conc[t, p] + model.v_dot_dil[t, p]
model.process_flow_balance = pyo.Constraint(model.T, model.processes, rule=process_flow_balance)



# limits inflow to each process
def inlet_flow_max(model, t, p):
    return model.v_dot_in[t, p] <= model.V_dot_max * model.y_process_active[p]
model.inlet_flow_max = pyo.Constraint(model.T, model.processes, rule=inlet_flow_max)
 


# limits seawater feed to one process 
def single_sw_feed(model):
    return sum(model.y_link_active['SW', p] for p in model.processes if ('SW', p) in model.links) <= 1
model.single_sw_feed = pyo.Constraint(rule=single_sw_feed)



# limits downstream process to one upstream process
def single_upstream(m, q):
    return sum(m.y_link_active[(p2,q)] for (p2,q2) in m.links if q2==q) <= 1
model.single_upstream = pyo.Constraint(model.processes, rule=single_upstream)



# check for downstream link and sends outlet flow to one downstream process
def single_downstream(m, p):
    terms = [m.y_link_active[(p,q)] for (p2,q) in m.links if p2 == p]
    if not terms:
        return pyo.Constraint.Skip
    return sum(terms) <= 1
model.single_downstream = pyo.Constraint(model.processes, rule=single_downstream)



# enforces link between processes can only be activated if the upstream process is activated
def upstream_activation(model, q, p):
    if q == 'SW':
        return pyo.Constraint.Skip
    return model.y_link_active[(q,p)] <= model.y_process_active[q]
model.upstream_activation = pyo.Constraint(model.links, rule=upstream_activation)



# enforces link between processes can only be activated if the downstream process is activated
def downstream_activation(model, q, p):
    return model.y_link_active[q,p] <= model.y_process_active[p]
model.downstream_activation = pyo.Constraint(model.links, rule=downstream_activation)



# relates inflow to process as sum of upstream links
def inlet_flow(m, t, p):
    inflow = sum(m.v_dot_link[t, (p2,p)] for (p2,q) in m.links if q==p)
    return m.v_dot_in[t, p] == inflow
model.inlet_flow = pyo.Constraint(model.T, model.processes, rule=inlet_flow)



# limits seawater intake to maximum
def seawater_intake(m, t):
    return sum(m.v_dot_link[t, ('SW',q)] for (p2,q) in m.links if p2=='SW') <= m.V_dot_max
model.seawater_intake = pyo.Constraint(model.T, rule=seawater_intake)



# sets ro recovery ratio relationship 
def recovery_ro(model, t):
    return model.v_dot_dil[t, 'RO'] == model.Recovery_ro * model.v_dot_in[t, 'RO']
model.recovery_ro = pyo.Constraint(model.T, rule=recovery_ro)



# enforces no ions in diluted ro stream
def ro_diluted(model, t, p, i):
    if p!= 'RO':
        return pyo.Constraint.Skip
    return model.f_ion_dil[t, 'RO', i] == 0.0
model.ro_diluted = pyo.Constraint(model.T, model.processes, model.ions, rule=ro_diluted)



# enforces correct recovery ratio for med 
def recovery_med(model, t):
    return model.v_dot_dil[t, 'MED'] == (
        model.v_dot_in[t, 'MED'] * model.Recovery_sw_to_med * (1 - model.y_link_active['RO','MED']) + 
        model.v_dot_in[t, 'MED'] * model.Recovery_ro_to_med * model.y_link_active['RO','MED'])
model.recovery_med = pyo.Constraint(model.T, rule=recovery_med)



# enforces no ions in diluted med stream
def med_diluted(model, t, p, i):
    if p!= 'MED':
        return pyo.Constraint.Skip
    return model.f_ion_dil[t, 'MED', i] == 0.0
model.med_diluted = pyo.Constraint(model.T, model.processes, model.ions, rule=med_diluted)



# enforces no concentrated flow from cry
def cry_no_concentrate(model, t):
    return model.v_dot_conc[t, 'CRY'] == 0.0
model.cry_no_concentrate = pyo.Constraint(model.T, rule=cry_no_concentrate)

    

# enforces no ions in diluted cry stream 
def cry_diluted(model, t, i):
    return model.f_ion_dil[t, 'CRY', i] == 0.0
model.cry_diluted = pyo.Constraint(model.T, model.ions, rule=cry_diluted)


# enforces ratio between na and li in ed stream
def ed_precip_ratio(model,t):
    return model.f_ion_dil[t,'ED','Na'] == model.Ratio_na_to_li * model.f_ion_dil[t,'ED','Li']
model.ed_precip_ratio = pyo.Constraint(model.T, rule=ed_precip_ratio)



# relates li recovered to yield 
def li_recovery(model,t):
    return model.li_recovered_ed[t] == model.yield_ed_li * model.f_ion_in[t,'ED','Li']
model.li_recovery = pyo.Constraint(model.T, rule=li_recovery)



# equates li recovered with amount in diluted stream
def li_mass_balance(model,t):
    return model.li_recovered_ed[t] == model.f_ion_dil[t,'ED','Li']
model.li_mass_balance = pyo.Constraint(model.T, rule=li_mass_balance)



# splits na between ed outlet streams
def ed_split(model,t):
    return (model.f_ion_conc[t,'ED','Na'] * model.v_dot_dil[t,'ED']
        == model.f_ion_dil[t,'ED','Na'] * model.v_dot_conc[t, 'ED'])
model.ed_split = pyo.Constraint(model.T, rule=ed_split)


# splits cl between ed outlet streams
def ed_cl_split(model,t):
    return (model.f_ion_conc[t,'ED','Cl'] * model.v_dot_dil[t,'ED']
            == model.f_ion_dil[t,'ED','Cl'] * model.v_dot_conc[t,'ED'])
model.ed_cl_split = pyo.Constraint(model.T, rule=ed_cl_split)



# sets capacity on processes
def capacity_limit_rule(model,t,p):
    return model.v_dot_in[t,p] <= model.v_dot_cap[p]/24
model.capacity_limit = pyo.Constraint(model.T, model.processes, rule=capacity_limit_rule)



# allows flow on link to be some multiple of seawater inflow
def ion_link_gate(model,t,p,q,i):
    M = pyo.value(model.Conc_sw[i]) * pyo.value(model.V_dot_max)
    return model.f_ion_link[t, (p,q), i] <= M * model.y_link_active[(p,q)]
model.ion_link_gate = pyo.Constraint(model.T, model.links, model.ions, rule=ion_link_gate)



############################ NEW CONSTRAINTS ################################

# limits link flow to less than upstream process concentrate flow:
def upper_link_flow_balance(model,t,p,q):
    if p == 'SW':
        return pyo.Constraint.Skip
    return model.v_dot_link[t,(p,q)] <= model.v_dot_conc[t,p]
model.upper_link_flow_balance = pyo.Constraint(model.T, model.links, rule=upper_link_flow_balance)


# requires link flow to be upstream process concentrate flow if active
def lower_link_flow_balance(model,t,p,q):
    if p == 'SW':
        return pyo.Constraint.Skip
    return model.v_dot_link[t,(p,q)] >= model.v_dot_conc[t,p] - model.V_dot_max*(1 - model.y_link_active[(p,q)])
model.lower_link_flow_balance = pyo.Constraint(model.T, model.links, rule=lower_link_flow_balance)



# requires incoming link for process to be active
def process_requires_incoming(model, p):
    if p == 'SW':
        return pyo.Constraint.Skip
    incoming = [u for (u,v) in model.links if v == p]
    if not incoming:
        return pyo.Constraint.Skip
    return model.y_process_active[p] <= sum(model.y_link_active[u,p] for u in incoming)
model.process_requires_incoming = pyo.Constraint(model.processes, rule=process_requires_incoming)



# creates big-M for each ion
M_ion       = {i:pyo.value(model.V_dot_max) * 20.0 * pyo.value(model.Conc_sw[i]) for i in model.ions}
model.M_ion = pyo.Param(model.ions, initialize=M_ion, within=pyo.NonNegativeReals)



# sets upper limit on ion concentrate along a link
def ion_link_upper(model,t,p,q,i):
    if p == 'SW': 
        return pyo.Constraint.Skip
    return model.f_ion_link[t,(p,q),i] <= model.f_ion_conc[t,p,i]
model.ion_link_upper  = pyo.Constraint(model.T, model.links, model.ions, rule=ion_link_upper)



# sets lower limit on ion concentration along a link
def ion_link_lower(model,t,p,q,i):
    if p == 'SW': 
        return pyo.Constraint.Skip
    return model.f_ion_link[t,(p,q),i] >= model.f_ion_conc[t,p,i] - model.M_ion[i] * (1 - model.y_link_active[(p,q)])
model.ion_link_lower  = pyo.Constraint(model.T, model.links, model.ions, rule=ion_link_lower)



# guarantees no ion flow on inactive links
def ion_link_off(model,t,p,q,i):
    if p == 'SW': 
        return pyo.Constraint.Skip
    return model.f_ion_link[t,(p,q),i] <= model.M_ion[i] * model.y_link_active[(p,q)]
model.ion_link_off = pyo.Constraint(model.T, model.links, model.ions, rule=ion_link_off)



V_dot_min_total = 1e-3
Price_salt      = 0.0

# enforces processes must be used if binary active
def process_used_if_active(model,p):
    if p == 'SW':
        return pyo.Constraint.Skip
    return sum(model.v_dot_in[t,p] for t in model.T) >= V_dot_min_total * model.y_process_active[p]
model.process_used_if_active = pyo.Constraint(model.processes, rule=process_used_if_active)



# enforces links must be used if binary active
def link_used_if_active(model,p,q):
    return sum(model.v_dot_link[t,(p,q)] for t in model.T) >= V_dot_min_total * model.y_link_active[(p,q)]
model.link_used_if_active = pyo.Constraint(model.links, rule=link_used_if_active)

 

# enforces downstream endpoint required for all links
def link_requires_downstream_rule(model,p,q):
    return model.y_link_active[p, q] <= model.y_process_active[q]
model.link_requires_downstream = pyo.Constraint(model.links, rule=link_requires_downstream_rule)



# enforces upstream point required for links
def link_requires_upstream_rule(model,p,q):
    if p == 'SW':
        return pyo.Constraint.Skip
    return model.y_link_active[p, q] <= model.y_process_active[p]
model.link_requires_upstream = pyo.Constraint(model.links, rule=link_requires_upstream_rule)


    
# enforces if process is on it must have at least one incoming link
def process_requires_incoming(m, p):
    if p == 'SW':
        return pyo.Constraint.Skip
    incoming = [q for (q,v) in m.links if v == p]
    if not incoming:
        return pyo.Constraint.Skip
    return m.y_process_active[p] <= sum(m.y_link_active[(q,p)] for q in incoming)

model.process_requires_incoming = pyo.Constraint(model.processes, rule=process_requires_incoming)



############################ NEW CONSTRAINTS ################################


annualize = 8760/len(model.T)

profit_yr = annualize * sum(
    model.Price_li * model.li_recovered_ed[t]
  + model.Price_water * (model.v_dot_dil[t,'RO'] + model.v_dot_dil[t,'MED'] + model.v_dot_dil[t,'CRY'])
  + model.Price_elec[t] * model.elec_sold[t] + Price_salt * sum(model.f_ion_in[t,'CRY', i] for i in model.ions)
  for t in model.T)



capex_total = sum(model.K_process[p] * model.v_dot_cap[p] for p in model.processes) + model.Cost_tes_cold*model.m_cold_tes_max + model.Cost_tes_hot*model.m_hot_tes_max 


capex_yr    = capex_total/30

model.obj   = pyo.Objective(expr=(profit_yr - capex_yr), sense=pyo.maximize)

                          



solver  = pyo.SolverFactory("gurobi")
solver.options['NonConvex'] = 2
results = solver.solve(model, tee=True)




model.solutions.load_from(results)

print("\n=== Process Status and Capacity ===")
for p in model.processes:
    y   = pyo.value(model.y_process_active[p])
    cap = pyo.value(model.v_dot_cap[p])
    print(f"{p:>4s} | ON = {int(round(y))} | Capacity = {cap:.3f}")
    
    
print("\n=== Active process links ===")
for (u, v) in model.links:
    if pyo.value(model.y_link_active[u, v]) > 0.5:
        print(f"{u:>4s} → {v:<4s}")
    

t = model.T.last()   # or any representative time

print(f"\n=== Flow & Ion Mass Snapshot at t = {t} ===")

for p in model.processes:
    if pyo.value(model.y_process_active[p]) < 0.5:
        continue

    print(f"\nProcess {p}")
    print(
        f"  Volumetric flows [m3/h]: "
        f"in = {pyo.value(model.v_dot_in[t,p]):8.2f} | "
        f"conc = {pyo.value(model.v_dot_conc[t,p]):8.2f} | "
        f"dil = {pyo.value(model.v_dot_dil[t,p]):8.2f}"
    )

    print("  Ion mass flows [kg/h]:")
    for i in model.ions:
        f_in   = pyo.value(model.f_ion_in[t,p,i])
        f_conc = pyo.value(model.f_ion_conc[t,p,i])
        f_dil  = pyo.value(model.f_ion_dil[t,p,i])

        print(
            f"    {i:>3s} | "
            f"in = {f_in:10.3e} | "
            f"conc = {f_conc:10.3e} | "
            f"dil = {f_dil:10.3e}")



print("\nProcesses with y=1 but zero capacity:")
for p in model.processes:
    if pyo.value(model.y_process_active[p]) > 0.5 and pyo.value(model.v_dot_cap[p]) < 1e-6:
        print(p)

print("\nLinks with y=1 but zero total flow:")
for (p,q) in model.links:
    if pyo.value(model.y_link_active[(p,q)]) > 0.5:
        tot = sum(pyo.value(model.v_dot_link[t2,(p,q)]) for t2 in model.T)
        if tot < 1e-6:
            print(p,"->",q)


# store data in csv 

out_dir = "results"
os.makedirs(out_dir, exist_ok=True)

import csv

with open("solution_summary.csv", "w", newline="") as f:
    writer = csv.writer(f)

    # -------------------------
    # PROCESS SUMMARY (static)
    # -------------------------
    writer.writerow(["# PROCESS SUMMARY"])
    writer.writerow(["process", "active", "capacity_m3_per_day"])

    for p in model.processes:
        writer.writerow([
            p,
            int(round(pyo.value(model.y_process_active[p]))),
            pyo.value(model.v_dot_cap[p])
        ])

    # blank line
    writer.writerow([])

    # -------------------------
    # TIME SERIES (dynamic)
    # -------------------------
    writer.writerow(["# TIME SERIES"])
    writer.writerow([
        "time",
        "m_dot_hpt_kg_s",
        "m_dot_extract_kg_s",
        "m_hot_tes_kg",
        "m_cold_tes_kg",
        "w_dot_gen_kWe"

    ])

    for t in model.T:
        writer.writerow([
            t,
            pyo.value(model.m_dot_hpt[t]),
            pyo.value(model.m_dot_extract[t]),
            pyo.value(model.m_hot_tes[t]),
            pyo.value(model.m_cold_tes[t]),
            pyo.value(model.w_dot_gen[t])
        ])
        
    writer.writerow([])
    writer.writerow(["# PROCESS FLOWS (m3/h)"])
    writer.writerow(["time", "process", "v_in", "v_conc", "v_dil"])
    
    for t in model.T:
        for p in model.processes:
            writer.writerow([
                t, p,
                pyo.value(model.v_dot_in[t, p]),
                pyo.value(model.v_dot_conc[t, p]),
                pyo.value(model.v_dot_dil[t, p]),
            ])


