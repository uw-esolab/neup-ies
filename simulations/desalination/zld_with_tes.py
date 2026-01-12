import pyomo.environ as pyomo
import pandas as pd
from pyomo.opt import SolverFactory


# calculations from med model 


N            = 8760                                          # units: --,      number of time steps (1 year)
Q_dot        = 950                                           # units: MWth,    reactor outpt
W_dot_nom    = 494                                           # units: MWe,     power cycle electricity generation
Effic_cycle  = 0.52                                          # units: --,      nominal efficiency of power cycle
W_dot_gen    = W_dot_nom                                     # units: MWe,     power generated in cycle, used in runs where turbine is oversized 


Cp_salt             = 1.5                                 # units: kJ/kg-K, specific heat of molten salt
T_salt_in           = 565                                 # units: deg C,   temperature of molten salt exiting reactor (to hot ht storage)
T_salt_out          = 330                                 # units: deg C,   temperature of molten salt after salt-to-steam hx (to cold ht storage)
Delta_T_salt        = T_salt_in - T_salt_out              # units: deg C,   change in molten salt temperature across salt-to-steam hx
M_dot_salt          = Q_dot / (Cp_salt * Delta_T_salt)    # units: kg/s,    molten salt flow from reactor


H_steam_cyc_in      = 1338                                # units: kJ/kg,   steam enthalpy entering salt-to-steam hx
H_steam_cyc_out     = 3464                                # units: kJ/kg,   steam enthalpy exiting salt-to-steam hx
Delta_H_steam_cyc   = H_steam_cyc_out - H_steam_cyc_in    # units: kJ/kg,   change in steam enthalpy across salt-to-steam hx
M_dot_steam_cyc_nom = Q_dot / Delta_H_steam_cyc           # units: kg/s,    nominal mass flow of steam across salt-to-steam hx


H_steam_ext_in      = 2747.9                              # units: kJ/kg,   steam enthalpy entering steam-to-water hx
H_steam_ext_out     = 209.7                               # units: kJ/kg,   steam enthalpy exiting steam-to-water hx
Delta_H_steam_ext   = H_steam_ext_in - H_steam_ext_out    # units: kJ/kg,   change in steam enthalpy across steam-to-water hx

Delta_H_water       = 421.3                               # units: kJ/kg,   change in water enthalpy across steam-to-water hx 
Frac_steam_ext_nom  = 0.53
M_dot_steam_ext_nom = M_dot_steam_cyc_nom * Frac_steam_ext_nom
M_dot_water_nom     = M_dot_steam_ext_nom * Delta_H_steam_ext / Delta_H_water # units: kg/s,   nominal mass flow of water across steam-to-water hx


Price_storage_salt       = 29.8   # units: $/kWh,      price of molten salt
Convert_kWh_to_kg_salt   = Price_storage_salt / 3600 * Cp_salt * Delta_T_salt


Price_storage_water      = 20.0   # units: $/kWh
Convert_kWh_to_kg_water  = Price_storage_water


Schedule_effic           = pd.read_csv('sched_eff.csv', usecols=['Efficiency']).to_numpy().flatten()

Schedule_elec            = pd.read_csv('DIABLOCN_2_N001_PRICES.csv', usecols=['LMP']).to_numpy().flatten()

Schedule_m_dot_salt      = [M_dot_salt for n in range(N)]



# create model
model = pyomo.ConcreteModel()



model.T = pyomo.Set(initialize=range(1, N+1))  # for example


# define sets
model.T         = pyomo.Set(initialize=range(1,N+1))
model.processes = pyomo.Set(initialize=['RO', 'MED', 'ED', 'CRY'])
model.links     = pyomo.Set(dimen=2, initialize=[
    ('SW', 'RO'),
    ('SW', 'MED'),
    ('SW', 'ED'),
    ('SW', 'CRY'),
    ('RO', 'MED'),
    ('RO', 'ED'),
    ('RO', 'CRY'),
    ('MED', 'ED'),
    ('MED', 'CRY'),
    ('ED', 'CRY')])


model.ions      = pyomo.Set(initialize=['Li', 'Na', 'Cl'])
#model.cycles    = pyomo.Set(initialize=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'])




# define parameters - med model

model.nt                      = pyomo.Param(initialize=N, domain=pyomo.Integers)                                # units: --,           number of time steps
model.Delta_t                 = pyomo.Param(initialize=1)                                                       # units: hr,           time increments
model.W_dot_max               = pyomo.Param(initialize=W_dot_nom)                                               # units: MWe,          maximum output from power cycle
model.W_dot_min               = pyomo.Param(initialize=0.25*W_dot_gen)                                          # units: MWe,          minimum output from power cycle, 25% of nominal value
model.W_dot_ramp_max          = pyomo.Param(initialize=0.1*W_dot_gen)                                           # units: MWe,          maximum change in power cycle output from one time step to the next
#model.Eff_cycle_nom           = pyomo.Param(initialize=0.52)                                                    # units: --,           power cycle efficiency
model.Eff_cycle               = pyomo.Param(model.T, initialize={t: float(Schedule_effic[t-1]) for t in model.T})           # units: --,           efficiency of power cycle as a function of ambient temperature
model.M_dot_salt_nom          = pyomo.Param(model.T, initialize = dict(zip(model.T, Schedule_m_dot_salt)))      # units: kg/s,         mass flow rate of salt through salt-to-steam hx
model.M_dot_steam_cyc_nom     = pyomo.Param(initialize=447.1)                                                   # units: kg/s,         nominal mass flow rate of steam through salt-to-steam hx
model.K_hx_salt_to_steam_cyc  = pyomo.Param(initialize=M_dot_steam_cyc_nom/M_dot_salt)                          # units: --,           steam flow through cycle per unit of salt flow 
model.K_hx_steam_ext_to_water = pyomo.Param(initialize=M_dot_water_nom/M_dot_steam_ext_nom)                     # units: --,           water flow through low-temperature tes per unit flow of extracted steam
model.K_hpt_power             = pyomo.Param(initialize=W_dot_gen/M_dot_steam_cyc_nom)                           # units: MW/kg/s,      power produced per unit flow of steam through the cycle (W_dot_nom/M_dot_steam_cyc_nom)
model.K_lpt_loss              = pyomo.Param(initialize=W_dot_gen/M_dot_steam_ext_nom)                           # units: MW/kg/s,      power lost per unit flow of extracted steam (W_dot_nom/M_dot_steam_ext_nom)
model.M_dot_hpt_max           = pyomo.Param(initialize=model.W_dot_max.value / model.K_hpt_power.value)
model.M_dot_hpt_min           = pyomo.Param(initialize=model.W_dot_min.value / model.K_hpt_power.value)
model.Frac_ext                = pyomo.Param(initialize=0.53)
model.Time_amor               = pyomo.Param(initialize=30)                                                      # units: yr,           time of amoritization                                                
model.Convert_time            = pyomo.Param(initialize=3600)                                                    # units: sec/hr        time converstion between seconds and hour
#model.Convert_med             = pyomo.Param(initialize=86.4)                                                    # units: m^3/day/kg/s  conversion between constant for med water production
model.M_lt_tes_init           = pyomo.Param(initialize=0)                                                       # units: kg,           initial storage inventory of low-temperature tes
model.M_ht_tes_init           = pyomo.Param(initialize=0)                                                       # units: kg,           initial storage inventory of high-temperature tes 
model.Cost_ramp_power         = pyomo.Param(initialize=43.75*W_dot_gen/M_dot_steam_cyc_nom)                     # units: $/ 
#model.Cost_ramp_med           = pyomo.Param(initialize=0.2)                     
model.Cost_ht_tes             = pyomo.Param(initialize=Price_storage_salt/model.Convert_time.value*Cp_salt*Delta_T_salt)     # units: $/kg      cost of high-temperature (molten salt) storage
model.Cost_lt_tes             = pyomo.Param(initialize=Price_storage_water/model.Convert_time.value*Delta_H_water)           # units: $/kg      cost of low-temperature (pressurized water) storage
model.Price_water             = pyomo.Param(initialize=0.0021)                                                  # units: $/kg         selling price of water 
model.Price_elec              = pyomo.Param(model.T, initialize={t: float(Schedule_elec[t-1]/1000) for t in model.T})            # units: $/kWh       selling price of electricity(divide by 1000 to get from MWh to kWh)
model.Slope_power_curve       = pyomo.Param(initialize=1.075)
#model.Slope_cost_med          = pyomo.Param(initialize=1095.3)
model.Interc_power_curve      = pyomo.Param(initialize=0.0748)
#model.Interc_cost_med         = pyomo.Param(initialize=2.05E07)
#model.Eff_med                 = pyomo.Param(initialize=0.70)                                                    # units: --,          water production per unit mass flow of heated water into med
model.F_lpt_post              = pyomo.Param(initialize=0.186)                                                   # units: --,          fraction of total LPT power generated downstream of extraction point


# define parameters - zld model

model.Recovery_ro             = pyomo.Param(initialize=0.50)   # units: --,     recovery ratio for ro fed seawater
model.Recovery_sw_to_med      = pyomo.Param(initialize=0.77)   # units: --,     recovery ratio for med fed seawater
model.Recovery_ro_to_med      = pyomo.Param(initialize=0.53)   # units: --,     recovery ratio for med fed ro concentrate
#model.Price_water             = pyomo.Param(initialize=1.0)    # units: $/m3,   sales price of desalinated water
#model.Price_elec              = pyomo.Param(initialize=0.1)    # units: $/kWh,  sales price of electricity
model.Price_li                = pyomo.Param(initialize=70.0) # units: $/kg,   sales price of lithium
model.V_dot_min                = pyomo.Param(initialize=1e-04)  # units: --,     little m used for logical constraints
model.Ratio_ed_na_li          = pyomo.Param(initialize=5.0)    # units: --,     optimal ratio between sodium and lithium in concentrated stream leaving ed for precipitation            
model.V_dot_max               = pyomo.Param(initialize=40000)  # units; m3/h,   maximum seawater inflow for desalination, based on current largest desalination plant in world 
#model.K_cycle                 = pyomo.Param(initialize=6000)   # units: $/kW,   capex for nuclear power plants
model.Conc_sw                 = pyomo.Param(model.ions, initialize={'Li': 0.00018, 'Na': 10.8, 'Cl': 19.3})                                                               # units: kg/m3     mass concentration of ions in seawater
model.Elec_required           = pyomo.Param(model.processes, initialize={'RO': 4.5, 'MED': 1.0, 'ED': 3.5,'CRY': 3.5})                                                # units: kWe/m3    electrical energy required per unit of feed 
model.Q_required              = pyomo.Param(model.processes, initialize={'RO': 0.0,  'MED': 25.0 ,'ED': 0.0, 'CRY': 90})                                                 # units: kW-th/m3  heat required per unit of feed (for ed and crystallization these values are calculated based on equations in the model)
model.K_process               = pyomo.Param(model.processes, initialize={'RO': 500.0,'MED': 900.0,'ED': 200.0, 'CRY': 950.0})                                                 # units: $/m3/day  capex for each process
#model.H_extract               = pyomo.Param(model.cycles, initialize={'C1': 2.82e3, 'C2': 2.79e3, 'C3': 2.49e3, 'C4': 2.38e3, 'C5': 2.32e3, 'C6': 2.21e3, 'C7': 2.20e3})  # units: kJ/kg     enthalpy of steam at each extraction point within the power cycle configurations - steam enthalpy going back into condenser
model.yield_ed_li             = pyomo.Param(initialize=1.0)    # units: --,   li recovery fraction, assumed 100%

# define variables - med model

model.m_ht_tes                = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg,            inventory in high-temp hot storage 
model.m_dot_ht_tes_out        = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg/s,          mass flow out of high-temp hot storage
model.m_dot_hpt               = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg/s,          mass flow at the high pressure turbine inlet 
model.m_dot_extract           = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg/s,          mass flow of steam extracted from power cycle
model.m_lt_tes                = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg,            inventory in low-temp hot storage
model.m_dot_lt_tes_out        = pyomo.Var(model.T, model.processes, domain=pyomo.NonNegativeReals)                       # units: kg/s,          mass flow out of low-temperature storage to med
model.w_dot_gen               = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: MW,            electricity produced by power cycle
#model.v_dot_med               = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg/s,          water produced by desal system
#model.v_dot_med_Delta_up      = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)                   # units: kg/s,          positive change in water produced relative to the previous time step
#model.v_dot_med_Delta_dn      = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)                   # units: kg/s,          negative change in water produced relative to the previous time step
model.m_ht_tes_max            = pyomo.Var(domain=pyomo.NonNegativeReals)                                # units: kWh,           high-temp thermal storage maximum inventory
model.m_lt_tes_max            = pyomo.Var(domain=pyomo.NonNegativeReals)                                # units: kWh,           low-temp thermal storage maximum inventory
#model.v_dot_med_max           = pyomo.Var(domain=pyomo.NonNegativeReals)                                # units: kg/s,          size of med system
#model.capex_med               = pyomo.Var(domain=pyomo.NonNegativeReals)                                # units: $              cost of med system
model.m_dot_hpt_Delta_up      = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg/s,          upward change in hpt mass flow from one time step to the next 
model.m_dot_hpt_Delta_dn      = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                       # units: kg/s,          downward change in hpt mass flow from one time step to the next
#model.y_med_on                = pyomo.Var(domain=pyomo.Binary)                                          # units: --,            binary activator if med system is turned on


# define variables - zld model

model.y_link_active           = pyomo.Var(model.links, within=pyomo.Binary)                              # units: --,    1 if there is flow from process q to process p; 0 otherwise
model.y_process_active        = pyomo.Var(model.processes, within=pyomo.Binary)                          # units: --,    1 if process p is active; 0 otherwise
#model.y_cycle_active          = pyomo.Var(model.cycles, within=pyomo.Binary)                             # units: --,    1 if cycle c is active; 0, otherwise
#model.y_q_source              = pyomo.Var(model.cycles, model.processes, within=pyomo.Binary)            # units: --,    1 if cycle c provides heat to process p; 0, otherwise
model.v_dot_in                = pyomo.Var(model.T, model.processes, within=pyomo.NonNegativeReals)                # units: m3/hr  volumetric flow rate entering each process
model.v_dot_conc              = pyomo.Var(model.T, model.processes, within=pyomo.NonNegativeReals)                # units: m3/hr  volumetric flow rate of concentrated stream leaving each process
model.v_dot_dil               = pyomo.Var(model.T, model.processes, within=pyomo.NonNegativeReals)                # units: m3/hr  volumetric flow rate of diluted stream leaving each process
#model.v_dot_sw_feed           = pyomo.Var(model.T, model.processes, within=pyomo.NonNegativeReals)                # units: m3/hr  volumetric flow rate of seawater entering each process 
model.v_dot_link              = pyomo.Var(model.T, model.links, within=pyomo.NonNegativeReals)                             # units: m3/hr  volumetric flow rate from an upstream process to a downstream process
#model.m_dot_cycle_unit        = pyomo.Var(model.cycles, bounds=(0, model.M_dot_max))                             # units: kg/s   steam flow rate through each unit cycle configuration   
#model.m_dot_extract_unit      = pyomo.Var(model.cycles, bounds=(0, 0.8*model.M_dot_max))                         # units: kg/s   steam extraction flow rate for each power cycle configuration
model.elec_used               = pyomo.Var(model.T, model.processes, within=pyomo.NonNegativeReals)                # units: kWh-e   electricity used by each process
#model.elec_generated_unit     = pyomo.Var(model.T, model.cycles, within=pyomo.NonNegativeReals)                   # units: kWh-e   electricity generated by each power cycle configuration
model.elec_sold               = pyomo.Var(model.T, within=pyomo.NonNegativeReals)                                 # units: kWh-e   electricity sold to the grid from each power cycle configuration   
#model.elec_allocated          = pyomo.Var(model.T, model.processes, within=pyomo.NonNegativeReals)                # units: kWh-e   electricity allocated from each power cycle to each process
#model.elec_generated_unmasked = pyomo.Var(model.T, model.cycles, within=pyomo.NonNegativeReals)
model.q_generated             = pyomo.Var(model.T, within=pyomo.NonNegativeReals)                   # units: kWh-th  heat generated by each power cycle from steam extractions
model.q_used                  = pyomo.Var(model.T, model.processes, within=pyomo.NonNegativeReals)                # units: kWh-th  heat used by each process
#model.q_allocated             = pyomo.Var(model.T, model.cycles, model.processes, within=pyomo.NonNegativeReals)  # units: kWh-th  heat allocated from each power cycle to each process
model.f_ion_in                = pyomo.Var(model.T, model.processes, model.ions, within=pyomo.NonNegativeReals)    # units: kg/h    flow of each ion into each process
#model.concentration_conc      = pyomo.Var(model.T, model.processes, model.ions, within=pyomo.NonNegativeReals)    # units: kg/h    mass concentration of each ion in the concentrated stream exiting each process
#model.concentration_dil       = pyomo.Var(model.T, model.processes, model.ions, within=pyomo.NonNegativeReals)    # units: kg/h    mass concentration of each ion in the diluted stream exiting each process
model.li_recovered_ed         = pyomo.Var(model.T, within=pyomo.NonNegativeReals)                                 # units: kg/h, total Li recovered


# new constraintsa

model.f_ion_conc = pyomo.Var(model.T, model.processes, model.ions, within=pyomo.NonNegativeReals)  # kg/h
model.f_ion_dil  = pyomo.Var(model.T, model.processes, model.ions, within=pyomo.NonNegativeReals)  # kg/h
model.f_ion_link = pyomo.Var(model.T, model.links,     model.ions, within=pyomo.NonNegativeReals)  # kg/h


def ion_conservation(m,t,p,i):
    if p == 'CRY':
        return pyomo.Constraint.Skip
    return m.f_ion_in[t,p,i] == m.f_ion_conc[t,p,i] + m.f_ion_dil[t,p,i]
model.ion_conservation = pyomo.Constraint(model.T, model.processes, model.ions, rule=ion_conservation)



def ion_mass_in(m,t,p,i):
    return m.f_ion_in[t,p,i] == sum(m.f_ion_link[t,(u,p),i]
                                   for (u,v) in m.links if v == p)
model.ion_mass_in = pyomo.Constraint(model.T, model.processes, model.ions, rule=ion_mass_in)


def sw_ion_link(m,t,p,i):
    if ('SW',p) not in m.links:
        return pyomo.Constraint.Skip
    return m.f_ion_link[t,('SW',p),i] == m.v_dot_link[t,('SW',p)] * m.Conc_sw[i]
model.sw_ion_link = pyomo.Constraint(model.T, model.processes, model.ions, rule=sw_ion_link)



# def sw_flow_gate(m,t,p):
#     if ('SW',p) not in m.links:
#         return pyomo.Constraint.Skip
#     return m.v_dot_link[t,('SW',p)] <= m.V_dot_max * m.y_link_active[('SW',p)]
# model.sw_flow_gate = pyomo.Constraint(model.T, model.processes, rule=sw_flow_gate)


# def link_gate(m,t,u,v):
#     if u == 'SW':
#         return pyomo.Constraint.Skip
#     return m.v_dot_link[t,(u,v)] <= m.V_dot_max * m.y_link_active[(u,v)]
# model.link_gate = pyomo.Constraint(model.T, model.links, rule=link_gate)


def link_gate_all(m, t, p, q):
    return m.v_dot_link[t,(p,q)] <= m.V_dot_max * m.y_link_active[(p,q)]

model.link_gate_all = pyomo.Constraint(model.T, model.links, rule=link_gate_all)


# skip enforcing flow balance if the process doesn't route concentrate anywhere
def conc_outflow_balance(m,t,source_proc):
    if source_proc == 'SW':
        return pyomo.Constraint.Skip
    outgoing_links = [
        (link_source, link_target)
        for (link_source, link_target) in m.links
        if link_source == source_proc and link_source != 'SW'
    ]
    if len(outgoing_links) == 0:
        return pyomo.Constraint.Skip
    return sum(m.v_dot_link[t, (source_proc, link_target)] for (_, link_target) in outgoing_links) == m.v_dot_conc[t, source_proc]

model.conc_outflow_balance = pyomo.Constraint(model.T, model.processes, rule=conc_outflow_balance)



# skip enforcing the rule if no links exist, otherwise enforce only one outgoing link 
def one_outgoing(m, t, source_proc):
    if source_proc == 'SW':
        return pyomo.Constraint.Skip
    outgoing_links = [
        (link_source, link_target)
        for (link_source, link_target) in m.links
        if link_source == source_proc and link_source != 'SW'
    ]
    if len(outgoing_links) == 0:
        return pyomo.Constraint.Skip
    return sum(m.y_link_active[(source_proc, link_target)] for (_, link_target) in outgoing_links) <= 1
model.one_outgoing = pyomo.Constraint(model.T, model.processes, rule=one_outgoing)



def route_conc_ions(m, t, source_proc, ion):
    if source_proc in ['SW', 'CRY']:
        return pyomo.Constraint.Skip
    outgoing_links = [
        (link_source, link_target)
        for (link_source, link_target) in m.links
        if link_source == source_proc and link_source != 'SW'
    ]
    if len(outgoing_links) == 0:
        return pyomo.Constraint.Skip
    return sum(m.f_ion_link[t, (source_proc, link_target), ion] for (_, link_target) in outgoing_links) == m.f_ion_conc[t, source_proc, ion]
model.route_conc_ions = pyomo.Constraint(model.T, model.processes, model.ions, rule=route_conc_ions)


# power cycle constraints

# convert steam to power in the power cycle
def constr_steam_to_power_conv(model, t):
    return model.w_dot_gen[t] == ((model.Eff_cycle[t] * W_dot_gen) / M_dot_steam_cyc_nom) * (model.Slope_power_curve * model.m_dot_hpt[t] - model.F_lpt_post * model.m_dot_extract[t] - model.Interc_power_curve * M_dot_steam_cyc_nom)
model.constr_steam_to_power_conv =  pyomo.Constraint(model.T, rule=constr_steam_to_power_conv)


# mass flow limit on extraction
def constr_extraction_limit(model, t):
    return model.m_dot_extract[t] <= model.Frac_ext * model.m_dot_hpt[t]
model.constr_extraction_limit = pyomo.Constraint(model.T, rule=constr_extraction_limit)


# limits flow to hpt based on size of turbine
def constr_bound_hpt_flow(model, t):
    return model.m_dot_hpt[t] <= model.M_dot_hpt_max
model.constr_bound_hpt_flow = pyomo.Constraint(model.T, rule=constr_bound_hpt_flow)



# relate salt mass flow out of cycle storage to steam mass flow in power cycle
def constr_salt_to_steam_conv(model,t):
    return model.m_dot_hpt[t] == model.K_hx_salt_to_steam_cyc * model.m_dot_ht_tes_out[t]
model.constr_salt_to_steam_conv = pyomo.Constraint(model.T, rule=constr_salt_to_steam_conv)



# maximum inventory limit on high-temp tes
def constr_ht_tes_max_inventory(model,t):
    return model.m_ht_tes[t] <= model.m_ht_tes_max
model.constr_ht_tes_max_inventory =  pyomo.Constraint(model.T, rule=constr_ht_tes_max_inventory)



# power ramping tracking upwards
def constr_hpt_track_up(model, t):
    return model.m_dot_hpt_Delta_up[t] >= model.m_dot_hpt[t] - model.m_dot_hpt[t-1]
model.constr_hpt_track_up = pyomo.Constraint(model.T-[1], rule=constr_hpt_track_up)



# power ramping tracking downwards
def constr_hpt_track_down(model, t):
    return model.m_dot_hpt_Delta_dn[t] >= model.m_dot_hpt[t-1] - model.m_dot_hpt[t]
model.constr_hpt_track_down = pyomo.Constraint(model.T-[1], rule=constr_hpt_track_down)



# power ramping limit upwards
def constr_cycle_ramp_up(model,t):
    return model.m_dot_hpt_Delta_up[t] <= model.W_dot_ramp_max / model.K_hpt_power * model.Delta_t
model.constr_cycle_ramp_up = pyomo.Constraint(model.T-[1], rule=constr_cycle_ramp_up)



# power ramping downwards
def constr_cycle_ramp_down(model,t):
    return model.m_dot_hpt_Delta_dn[t] <= model.W_dot_ramp_max / model.K_hpt_power * model.Delta_t
model.constr_cycle_ramp_down = pyomo.Constraint(model.T-[1], rule=constr_cycle_ramp_down)  



# mass conservation in high-temperature tes
def constr_ht_tes_inventory(model, t):
    
    if t > 1:
        return model.m_ht_tes[t] == (model.M_dot_salt_nom[t] - model.m_dot_ht_tes_out[t]) * model.Convert_time * model.Delta_t + model.m_ht_tes[t-1]
    
    else:
        return model.m_ht_tes[t] == (model.M_dot_salt_nom[t] - model.m_dot_ht_tes_out[t]) * model.Convert_time * model.Delta_t + model.M_ht_tes_init

model.constr_ht_tes_inventory = pyomo.Constraint(model.T, rule=constr_ht_tes_inventory)



# electricity used by each process is related to the electrical requirements of each process
def electricity_requirement_rule(model, t, p):
    return model.elec_used[t,p] == model.Elec_required[p] * model.v_dot_in[t,p]
model.electricity_requirement = pyomo.Constraint(model.T, model.processes, rule=electricity_requirement_rule)



# electricity generated must either be sold or allocated to a process
def electricity_balance_rule(model, t):
    return model.elec_sold[t] == model.w_dot_gen[t] - sum(model.elec_used[t, p] for p in model.processes)
model.electricity_balance = pyomo.Constraint(model.T, rule=electricity_balance_rule)



# discharge from lt-storage can either go to med or cry
def constr_lt_tes_heat_balance(model, t, p):
    return model.m_dot_lt_tes_out[t, p] * Delta_H_water * model.Delta_t == model.q_used[t, p]
model.constr_lt_tes_heat_balance = pyomo.Constraint(model.T, model.processes, rule=constr_lt_tes_heat_balance)



# heat used by each process must equal the kWh-th required
def constr_q_used_defn(model, t, p):
    return model.q_used[t, p] == model.Q_required[p] * model.v_dot_in[t, p]
model.constr_q_used_defn = pyomo.Constraint(model.T, model.processes, rule=constr_q_used_defn)



# desalination constraints

# conservation of volumetric flow through the inlet and outlet of each process
def process_flow_balance_rule(model, t, p):
    return model.v_dot_in[t, p] == model.v_dot_conc[t, p] + model.v_dot_dil[t, p]
model.process_flow_balance = pyomo.Constraint(model.T, model.processes, rule=process_flow_balance_rule)



# if a link is turned on there has to be flow on it
def logical_no_flow_means_inactive_rule(model, t, q, p):
    return model.v_dot_link[t, q, p] >= model.V_dot_min * model.y_link_active[q, p]
model.logical_no_flow_means_inactive = pyomo.Constraint(model.T, model.links, rule=logical_no_flow_means_inactive_rule)



# process can only have inlet flow if its binary activation is turned on
def logical_inflow_maximum_activation_rule(model, t, p):
    return model.v_dot_in[t, p] <= model.V_dot_max * model.y_process_active[p]
model.logical_inflow_maximum_activation = pyomo.Constraint(model.T, model.processes, rule=logical_inflow_maximum_activation_rule)



# process must have inlet flow if it's binary activation is turned on
def logical_inflow_minimum_activation_rule(model, t, p):
    return model.v_dot_in[t, p] >= model.V_dot_min * model.y_process_active[p]
model.logical_inflow_minimum_activation = pyomo.Constraint(model.T, model.processes, rule=logical_inflow_minimum_activation_rule)
   


# seawater feed can only go to one process to start the chain of processes
def logical_single_feed_target_rule(model):
    return sum(model.y_link_active['SW', p] for p in model.processes if ('SW', p) in model.links) <= 1
model.logical_single_feed_target = pyomo.Constraint(rule=logical_single_feed_target_rule)



# each downstream process can be fed by only one upstream process
def logical_single_upstream_rule(model, p):
    return sum(model.y_link_active[upstream, p] for (upstream, downstream) in model.links if downstream==p) <= 1
model.logical_single_upstream = pyomo.Constraint(model.processes, rule=logical_single_upstream_rule)



# each process can send outlet flow to only one downstream process, checks to make sure there is a downstream process and if not skips the constraint
def logical_single_downstream_rule(model, p):
    terms = [model.y_link_active[p, downstream] for (upstream, downstream) in model.links if upstream == p]
    if not terms:
        return pyomo.Constraint.Skip
    return sum(terms) <= 1
model.logical_single_downstream = pyomo.Constraint(model.processes, rule=logical_single_downstream_rule)



# link between processes can only be activated if the upstream process is activated
def logical_link_activation_upstream_rule(model, upstream, downstream):
    if upstream == 'SW':
        return pyomo.Constraint.Skip
    return model.y_link_active[upstream, downstream] <= model.y_process_active[upstream]
model.logical_link_activation_upstream = pyomo.Constraint(model.links, rule=logical_link_activation_upstream_rule)



# link between processes can only be activated if the downstream process is activated
def logical_link_activation_downstream_rule(model, upstream, downstream):
    return model.y_link_active[upstream, downstream] <= model.y_process_active[downstream]
model.logical_link_activation_downstream = pyomo.Constraint(model.links, rule=logical_link_activation_downstream_rule)



# the flow from one link to another is 0 if the link is inactive; otherwise it is limited by a big M
def logical_link_capacity_rule(model, t, upstream, downstream):
    return model.v_dot_link[t, upstream, downstream] <= model.V_dot_max * model.y_link_active[upstream, downstream]
model.logical_link_capacity = pyomo.Constraint(model.T, model.links, rule=logical_link_capacity_rule)



# # the ion flow to each process is related to what is feeding it 
# def logical_ion_mass_in_rule(model, t, p, i):
    
#     seawater_term = 0.0
#     if ('SW', p) in model.links:
#         seawater_term = model.v_dot_sw_feed[t, p] * model.Conc_sw[i] * model.y_link_active['SW', p]
           
#     upstream_term = sum(
#         model.v_dot_link[t, upstream, p] * model.concentration_conc[t, upstream, i]
#         for (upstream, downstream) in model.links if downstream == p and upstream != 'SW'
#     )
    
#     return model.f_ion_in[t, p, i] == seawater_term + upstream_term
# model.logical_ion_mass_in = pyomo.Constraint(model. T, model.processes, model.ions, rule=logical_ion_mass_in_rule)




# # ensures that the flow on the links is correct depending on the processes being used
# def logical_link_flow_balance(model, t, upstream, p):
    
#     if (upstream, p) not in model.links or upstream == 'SW':
#         return pyomo.Constraint.Skip

#     else:
#         return model.v_dot_link[t, upstream, p] == model.v_dot_conc[t, upstream] * model.y_link_active[upstream, p]
    
# model.logical_link_flow_balance = pyomo.Constraint(model.T, model.links, rule=logical_link_flow_balance)



# ion conservation from inlet to concentrated and diluted stream
# def ion_conservation_rule(model, t, p, i):
#     if p == 'CRY':
#         return pyomo.Constraint.Skip
#     else:
#         return model.f_ion_in[t, p, i] == model.v_dot_conc[t, p] * model.concentration_conc[t, p, i] + model.v_dot_dil[t, p] * model.concentration_dil[t, p, i]
# model.ion_conservation = pyomo.Constraint(model.T, model.processes, model.ions, rule=ion_conservation_rule)



# the volumetric flow into a process is the sum of the flow on the links connected to the process upstream
def inlet_flow_balance_rule(model, t, p):
    inflow_from_links = sum(model.v_dot_link[t, q, downstream] for (q, downstream) in model.links if downstream == p)
    return model.v_dot_in[t, p] == inflow_from_links
model.inlet_flow_balance = pyomo.Constraint(model.T, model.processes, rule=inlet_flow_balance_rule)



# seawater intake flow is limited by parameter, sums over all connections that include sw as a source node, however a single upstream and downstream node is also enforced in the model
def feed_intake_capacity_rule(model, t):
    return sum(model.v_dot_link[t, 'SW', p] for (upstream, p) in model.links if upstream == 'SW') <= model.V_dot_max
model.feed_intake_capacity = pyomo.Constraint(model.T, rule=feed_intake_capacity_rule)



# diluted (potable) water is a function of recovery ratio
def recovery_ratio_ro_rule(model, t):
    return model.v_dot_dil[t, 'RO'] == model.Recovery_ro * model.v_dot_in[t, 'RO']
model.recovery_ratio_ro = pyomo.Constraint(model.T, rule=recovery_ratio_ro_rule)


# no ions in diluted ro stream
def diluate_concentration_ro_rule(model, t, p, i):
    if p!= 'RO':
        return pyomo.Constraint.Skip
    return model.f_ion_dil[t, 'RO', i] == 0.0
model.diluate_concentration_ro = pyomo.Constraint(model.T, model.processes, model.ions, rule=diluate_concentration_ro_rule)



# diluted water is a function of recovery ratio - recovery ratio for med depends on if it is fed by seawater or ro brine
def logical_recovery_ratio_med_rule(model, t):
    return model.v_dot_dil[t, 'MED'] == (
        model.v_dot_in[t, 'MED'] * model.Recovery_sw_to_med * (1 - model.y_link_active['RO','MED']) + 
        model.v_dot_in[t, 'MED'] * model.Recovery_ro_to_med * model.y_link_active['RO','MED'])
    model.logical_recovery_ratio_med = pyomo.Constraint(model.T, rule=logical_recovery_ratio_med_rule)


# no ions in diluted med stream
def diluate_concentration_med_rule(model, t, p, i):
    if p!= 'MED':
        return pyomo.Constraint.Skip
    return model.f_ion_dil[t, 'MED', i] == 0.0
model.diluate_concentration_med = pyomo.Constraint(model.T, model.processes, model.ions, rule=diluate_concentration_med_rule)



# low-temp tes inventory tracking between time steps
def constr_lt_tes_track(model, t):
    if t > 1:
        return model.m_lt_tes[t] == model.m_lt_tes[t-1] + (
            model.K_hx_steam_ext_to_water * model.m_dot_extract[t] 
            - sum(model.m_dot_lt_tes_out[t, p] for p in ['MED', 'CRY'])
        ) * model.Delta_t
    else:
        return model.m_lt_tes[t] == model.M_lt_tes_init + (
            model.K_hx_steam_ext_to_water * model.m_dot_extract[t] 
            - sum(model.m_dot_lt_tes_out[t, p] for p in ['MED', 'CRY'])
        ) * model.Delta_t

model.constr_lt_tes_track = pyomo.Constraint(model.T, rule=constr_lt_tes_track)



# I think this one will no longer be true since it's in ion mass terms not concentration terms chloride ions are split between the two streams
# def ed_cl_concentration_rule(model, t):
#     return model.f_ion_conc[t, 'ED', 'Cl'] == model.f_ion_dil[t, 'ED', 'Cl']
# model.ed_cl_concentration = pyomo.Constraint(model.T, rule=ed_cl_concentration_rule)


# crystallization constraints

# there is no concentrated brine from crystallization, only potable water and solid salts
def cry_no_concentrated_flow_rule(model, t):
    return model.v_dot_conc[t, 'CRY'] == 0.0
model.cry_no_concentrated_flow = pyomo.Constraint(model.T, rule=cry_no_concentrated_flow_rule)

    
# there are no ions in the diluted stream outflowing from crystallization
def diluate_concentration_cry_rule(model, t, i):
    return model.f_ion_dil[t, 'CRY', i] == 0.0
model.diluate_concentration_cry = pyomo.Constraint(model.T, model.ions, rule=diluate_concentration_cry_rule)



# if a process is inactive, the concentration of ions in the concentrated stream must equal 0 (added because some ions were being recorded in processes that were turned off)
def concentration_conc_inactive_rule(model, t, p, i):
    return model.f_ion_conc[t,p,i] <= 1e9 * model.y_process_active[p]
model.concentration_conc_inactive = pyomo.Constraint(model.T, model.processes, model.ions, rule=concentration_conc_inactive_rule)


# if a process is inactive, the concentration of ions in the diluted stream must equal 0
def concentration_dil_inactive_rule(model, t, p, i):
    return model.f_ion_dil[t,p,i] <= 1e9 * model.y_process_active[p]
model.concentration_dil_inactive = pyomo.Constraint(model.T, model.processes, model.ions, rule=concentration_dil_inactive_rule)


# links seawater feed to flow on the link between seawater and process p
# def sw_feed_link_coupling_rule(m, p):
#     if ('SW', p) not in m.links:
#         return pyo.Constraint.Skip
#     return m.v_dot_sw_feed[p] == m.v_dot_link['SW', p]
# model.sw_feed_link_coupling = pyo.Constraint(model.processes, rule=sw_feed_link_coupling_rule)


# solves for concentration of feed into ED 
#BL - when I remove this constraint I get weird behavior. Dont know why but leaving it in
# replace with constraint below
# def ed_feed_concentration_rule(m, i):
#     return m.f_ion_in['ED', i] == m.v_dot_in['ED'] * m.conc_feed['ED', i]
# model.ed_feed_concentration = pyo.Constraint(model.ions, rule=ed_feed_concentration_rule)

# I don't think this constraint is necessary
# def ed_ion_inflow_rule(model, t, i):
#     return model.f_ion_in[t, 'ED', i] == sum(
#         model.f_ion_link[t, q, 'ED', i] for (q, p) in model.links if p == 'ED'
#     )
# model.ed_ion_inflow = pyomo.Constraint(model.T, model.ions, rule=ed_ion_inflow_rule)



 # lithium precipitation criteria is that concentrated ed stream needs 5 times as much sodium as lithium
# this has been changed from "equal" to "less than or equal to"
# changed to be on mass basis
def ed_precipitation_ratio_rule(model,t):
    return model.f_ion_dil[t,'ED','Na'] == model.Ratio_ed_na_li * model.f_ion_dil[t,'ED','Li']
model.ed_precipitation_ratio = pyomo.Constraint(model.T, rule=ed_precipitation_ratio_rule)


# lithium recovered is yield * what is entering ED
def li_recovery_rule(model,t):
    return model.li_recovered_ed[t] == model.yield_ed_li * model.f_ion_in[t,'ED','Li']
model.li_recovery = pyomo.Constraint(model.T, rule=li_recovery_rule)


# connects lithium recovered to diluted stream
# changed to be on mass balance
# def ed_li_conc_mass_balance_rule(m):
#     return m.li_recovered_ed == m.concentration_dil['ED','Li'] * m.v_dot_dil['ED']
# model.ed_li_conc_mass_balance = pyo.Constraint(rule=ed_li_conc_mass_balance_rule)


def ed_li_mass_balance_rule(model, t):
    return model.li_recovered_ed[t] == model.f_ion_dil[t, 'ED', 'Li']
model.ed_li_mass_balance = pyomo.Constraint(model.T, rule=ed_li_mass_balance_rule)



# connects lithium recovered to what is still in diluted stream
# changed to be on mass basis below
# def ed_li_dil_mass_balance_rule(m):
#     return m.concentration_conc['ED','Li'] * m.v_dot_conc['ED'] == m.f_ion_in['ED','Li'] - m.li_recovered_ed
# model.ed_li_dil_mass_balance = pyo.Constraint(rule=ed_li_dil_mass_balance_rule)

## not needed because already defined elsewhere
# def ed_li_mass_balance_rule(model, t):
#     return model.f_ion_conc[t, 'ED', 'Li'] == model.f_ion_in[t, 'ED', 'Li'] - model.li_recovered_ed[t]
# model.ed_li_mass_balance = pyomo.Constraint(model.T, rule=ed_li_mass_balance_rule)



# this makes the concentration in both sides of the ED stream equal...not sure if this is the best choice 
# (although the concentrations would be equal entering the stream,
# but without setting some type of relationship it will make one side extremely concentrated)
#BL thinks that while this constraint does nto affect the results, it should be left in to keep Na stream tidy.
# this is edited below
# def ed_na_conc_mass_balance_rule(m):
#     return m.concentration_conc['ED','Na'] == m.concentration_dil['ED', 'Na']
# model.ed_na_conc_mass_balance = pyo.Constraint(rule=ed_na_conc_mass_balance_rule)


# # this is still bilinear...
def ed_na_mass_proportional_split_rule(model, t):
    return (
        model.f_ion_conc[t, 'ED', 'Na'] * model.v_dot_dil[t, 'ED']
        ==
        model.f_ion_dil[t, 'ED', 'Na'] * model.v_dot_conc[t, 'ED']
    )
model.ed_na_mass_proportional_split = pyomo.Constraint(model.T, rule=ed_na_mass_proportional_split_rule)





def obj_rule(model):
    return sum(model.Price_elec[t] * model.elec_sold[t] for t in model.T)

model.obj = pyomo.Objective(rule=obj_rule, sense=pyomo.maximize)


solver = SolverFactory('gurobi') 

results = solver.solve(model, tee=True)


for t in model.T:
    print(t,
          pyomo.value(Schedule_elec[t]),
          pyomo.value(model.w_dot_gen[t]),
          pyomo.value(model.elec_sold[t]),
          pyomo.value(model.m_dot_hpt[t]),
          pyomo.value(model.m_dot_extract[t]))


print(results.solver.termination_condition)





