# import necessary packages
import pyomo.environ as pyo

"""Important Notes
The Crystallizer function was seemingly causing local minima and could result in garbage results. More details below 
but TLDR - use a constant crystallizer heat requirement.
Model is still nonlinear but behavior is now OK
"""

"""TEST CASES

Note: I could not resolve the base case heat/electricity estimates so picked new ones from (frankly) ChatGPT. This obviously
needs revisiting. However, the base estimates can be used to validate model behavior

Base case: SW->MED. Note: extremely good (but not implausible) heat diversion from power conversion cycle is driving this being best
                    Note2: not fully MED, it maxes out at a cycle extraction tradeoff point, which is interesting and good

All perturbations return to nominal before making new change

Note: at the 

Perturbation 1: Drop RO energy use to 3.0 - RO enters the system, MED is still present

Perturbation 2: 10x MED electricity use - RO only at maximum

Perturbation 3: water is worthless - no active links

Perturbation 4: lithium price x100. SW->MED->ED. Behavior is sensible

Perturbation 5: 10x MED electricity usage AND lithium price x100. SW->RO->MED->ED. MED is still used to help concentrate for ED

Perturbation 6: 100x MED electricity usage AND lithium price x100. SW->RO->ED. Now MED is not worth it

Perturbation 7: Salt is now worth $100/kg. Processed food becomes far too expensive. SW->RO->MED->CRY. No point in ED as Li is too cheap and it doesnt concentrate. But use CRY

Perturbation 8: Salt is worth $100/kg, Li price x100. Seawater mining takes off. RDO-3 team wins Nobel Prize for Engineering. SW->RO->MED->ED->CRY is built. $724Bn in revenue.

Perturbation 9: 10x MED electricity use, RO up to 4.6 electricity. No active links. This is important cos RO is v.marginal at baseline

Final note: there are still some differnces in the ED rules between 
"""


# build model

model = pyo.ConcreteModel()


# SETS 

# set of processes available for use
model.processes = pyo.Set(initialize=['RO', 'MED', 'ED', 'CRY'])


# set of all possible links between processes
model.links     = pyo.Set(dimen=2, initialize=[
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


# set of all ions considered
model.ions      = pyo.Set(initialize=['Li', 'Na', 'Cl'])


# set of all power cycle configurations
model.cycles    = pyo.Set(initialize=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'])


## PARAMETERS

model.Recovery_ro         = pyo.Param(initialize=0.50)   # units: --,     recovery ratio for ro fed seawater
model.Recovery_sw_to_med  = pyo.Param(initialize=0.77)   # units: --,     recovery ratio for med fed seawater
model.Recovery_ro_to_med  = pyo.Param(initialize=0.53)   # units: --,     recovery ratio for med fed ro concentrate
model.M_dot_max           = pyo.Param(initialize=67.07)  # units: kg/s,   maximum mass flow of steam through power cycle
model.Rho                 = pyo.Param(initialize=1025)   # units: kg/m3,  average density of seawater
model.Lambda              = pyo.Param(initialize=2257)   # units: kJ/kg,  latent heat of water at 100 degrees C
model.Time_conversion     = pyo.Param(initialize=3600)   # units: sec/hr, conversion between seconds and hours
model.Price_water         = pyo.Param(initialize=1.0)    # units: $/m3,   sales price of desalinated water
model.Price_elec          = pyo.Param(initialize=0.1)    # units: $/kWh,  sales price of electricity
model.Price_li            = pyo.Param(initialize=70.0) # units: $/kg,   sales price of lithium
model.Price_salt          = pyo.Param(initialize=0.0)    # units: $/kg    sales price of general salts
model.X_max               = pyo.Param(initialize=1)      # units: --,     maximum allowable number of units built for each power cycle configuration
model.Z_ed                = pyo.Param(initialize=0.9)    # units: --,     fraction of water evaporated from ED concentrate stream
model.Big_m               = pyo.Param(initialize=40000)  # units: --,     big M used for logical constraints
model.Big_m_heat          = pyo.Param(initialize=2e05)   # units: --,     big M used for logical heating constraints
model.Little_m            = pyo.Param(initialize=1e-05)  # units: --,     little m used for logical constraints
model.Ratio_ed_na_li      = pyo.Param(initialize=5.0)    # units: --,     optimal ratio between sodium and lithium in concentrated stream leaving ed for precipitation            
model.V_dot_max           = pyo.Param(initialize=40000)  # units; m3/h,   maximum seawater inflow for desalination, based on current largest desalination plant in world 
model.K_cycle             = pyo.Param(initialize=6000)   # units: $/kW,   capex for nuclear power plants
model.T_hour              = pyo.Param(initialize=1)      # units: hr,     need to convert from kW to kWh

#BL are the heat and electricity consumptions in kWh/m3 ???? (electric or thermal)

model.Conc_sw             = pyo.Param(model.ions, initialize={'Li': 0.00018, 'Na': 10.8, 'Cl': 19.3})                                                               # units: kg/m3     mass concentration of ions in seawater
model.Elec_required       = pyo.Param(model.processes, initialize={'RO': 4.5, 'MED': 1.0, 'ED': 3.5,'CRY': 3.5})                                                # units: kWe/m3    electrical energy required per unit of feed 
model.Q_required          = pyo.Param(model.processes, initialize={'RO': 0.0,  'MED': 25.0 ,'ED': 0.0, 'CRY': 90})                                                 # units: kW-th/m3  heat required per unit of feed (for ed and crystallization these values are calculated based on equations in the model)
model.K_process           = pyo.Param(model.processes, initialize={'RO': 500.0,'MED': 900.0,'ED': 200.0, 'CRY': 950.0})                                                 # units: $/m3/day  capex for each process
model.H_extract           = pyo.Param(model.cycles, initialize={'C1': 2.82e3, 'C2': 2.79e3, 'C3': 2.49e3, 'C4': 2.38e3, 'C5': 2.32e3, 'C6': 2.21e3, 'C7': 2.20e3})  # units: kJ/kg     enthalpy of steam at each extraction point within the power cycle configurations - steam enthalpy going back into condenser


# VARIABLES

model.y_link_active           = pyo.Var(model.links, within=pyo.Binary)                              # units: --,    1 if there is flow from process q to process p; 0 otherwise
model.y_process_active        = pyo.Var(model.processes, within=pyo.Binary)                          # units: --,    1 if process p is active; 0 otherwise
model.y_cycle_active          = pyo.Var(model.cycles, within=pyo.Binary)                             # units: --,    1 if cycle c is active; 0, otherwise
model.y_q_source              = pyo.Var(model.cycles, model.processes, within=pyo.Binary)            # units: --,    1 if cycle c provides heat to process p; 0, otherwise


model.v_dot_in                = pyo.Var(model.processes, within=pyo.NonNegativeReals)                # units: m3/hr  volumetric flow rate entering each process
model.v_dot_conc              = pyo.Var(model.processes, within=pyo.NonNegativeReals)                # units: m3/hr  volumetric flow rate of concentrated stream leaving each process
model.v_dot_dil               = pyo.Var(model.processes, within=pyo.NonNegativeReals)                # units: m3/hr  volumetric flow rate of diluted stream leaving each process
model.v_dot_sw_feed           = pyo.Var(model.processes, within=pyo.NonNegativeReals)                # units: m3/hr  volumetric flow rate of seawater entering each process 
model.v_dot_link              = pyo.Var(model.links, within=pyo.NonNegativeReals)                    # units: m3/hr  volumetric flow rate from an upstream process to a downstream process


model.m_dot_cycle_unit        = pyo.Var(model.cycles, bounds=(0, model.M_dot_max))                   # units: kg/s   steam flow rate through each unit cycle configuration   
model.m_dot_extract_unit      = pyo.Var(model.cycles, bounds=(0, 0.8*model.M_dot_max))               # units: kg/s   steam extraction flow rate for each power cycle configuration


model.elec_used               = pyo.Var(model.processes, within=pyo.NonNegativeReals)                # units: kWh-e   electricity used by each process
model.elec_generated_unit     = pyo.Var(model.cycles, within=pyo.NonNegativeReals)                   # units: kWh-e   electricity generated by each power cycle configuration
model.elec_sold               = pyo.Var(model.cycles, within=pyo.NonNegativeReals)                   # units: kWh-e   electricity sold to the grid from each power cycle configuration   
model.elec_allocated          = pyo.Var(model.cycles, model.processes, within=pyo.NonNegativeReals)  # units: kWh-e   electricity allocated from each power cycle to each process
model.elec_generated_unmasked = pyo.Var(model.cycles, within=pyo.NonNegativeReals)

model.q_generated             = pyo.Var(model.cycles, within=pyo.NonNegativeReals)                   # units: kWh-th  heat generated by each power cycle from steam extractions
model.q_used                  = pyo.Var(model.processes, within=pyo.NonNegativeReals)                # units: kWh-th  heat used by each process
model.q_allocated             = pyo.Var(model.cycles, model.processes, within=pyo.NonNegativeReals)  # units: kWh-th  heat allocated from each power cycle to each process


model.f_ion_in                = pyo.Var(model.processes, model.ions, within=pyo.NonNegativeReals)    # units: kg/h    flow of each ion into each process
model.concentration_conc      = pyo.Var(model.processes, model.ions, within=pyo.NonNegativeReals)    # units: kg/h    mass concentration of each ion in the concentrated stream exiting each process
model.concentration_dil       = pyo.Var(model.processes, model.ions, within=pyo.NonNegativeReals)    # units: kg/h    mass concentration of each ion in the diluted stream exiting each process



# PIECEWISE FUNCTIONS


# CONSTRAINTS

## NEW ADDITIONS TO MODEL


# new variables


model.conc_feed = pyo.Var(model.processes, model.ions, within=pyo.NonNegativeReals)
model.conc_feed['ED','Li'].setlb(0)
model.conc_feed['ED','Li'].setub(0.002) 

model.yield_ed_li     = pyo.Param(initialize=1.0)                  # units: --,   li recovery fraction, assumed 100%. 
model.li_recovered_ed = pyo.Var(within=pyo.NonNegativeReals)    # units: kg/h, total Li recovered


######## 


model.cry_salinity    = pyo.Var(domain=pyo.NonNegativeReals)
model.q_specific_cry  = pyo.Var(domain=pyo.NonNegativeReals)

model.power_slope     = pyo.Param(model.cycles, initialize= {'C1': -873, 'C2': -829, 'C3': -540, 'C4': -427, 'C5': -335, 'C6': -229, 'C7': -211})
model.power_intercept = pyo.Param(initialize=52920)

def cry_salinity_definition_rule(m):
    
    total_ion_flow = sum(m.f_ion_in['CRY', i] for i in m.ions)
    return m.cry_salinity * m.v_dot_in['CRY'] == total_ion_flow

model.cry_salinity_definition = pyo.Constraint(rule=cry_salinity_definition_rule)


"""BL: I believe this constraints are good but are messing up the solver. Instead have fixed this to 15 which 
I believe roughly lines up with the exit salinity of MED (or for that matter ED). I believe this is OK
because RO->CRY will never make sense, and all we are doing is slightly penalizing that
I cannot resolve the numbers from the paper cited with that for MED and have instead used a rule of thumb from online
that crystallizer needs about 4x as much heat and 3x as much power. This needs revisiting


def cry_specific_energy_rule(m):
    
    A = -0.0006
    B = 0.2648
    C = 2.2488

    return m.q_specific_cry ==  C  + (A * m.cry_salinity * m.cry_salinity) + (B * m.cry_salinity) 

model.cry_specific_energy = pyo.Constraint(rule=cry_specific_energy_rule)


def cry_total_energy_rule(m):    
    return m.q_used['CRY'] == (m.q_specific_cry * m.Rho / m.Time_conversion) * m.v_dot_in['CRY']
model.cry_total_energy = pyo.Constraint(rule=cry_total_energy_rule)
"""


def linear_power_generation_rule(m,c):
    return m.elec_generated_unmasked[c] == m.power_slope[c] * m.m_dot_extract_unit[c]  + m.power_intercept
model.linear_power_generation = pyo.Constraint(model.cycles, rule=linear_power_generation_rule)


########


# new constraints

# if a process is inactive, the concentration of ions in the concentrated stream must equal 0 (added because some ions were being recorded in processes that were turned off)
def concentration_conc_inactive_rule(m, p, i):
    return m.concentration_conc[p,i] <= 1e9 * m.y_process_active[p]
model.concentration_conc_inactive = pyo.Constraint(model.processes, model.ions, rule=concentration_conc_inactive_rule)


# if a process is inactive, the concentration of ions in the diluted stream must equal 0
def concentration_dil_inactive_rule(m, p, i):
    return m.concentration_dil[p,i] <= 1e9 * m.y_process_active[p]
model.concentration_dil_inactive = pyo.Constraint(model.processes, model.ions, rule=concentration_dil_inactive_rule)


# links seawater feed to flow on the link between seawater and process p
def sw_feed_link_coupling_rule(m, p):
    if ('SW', p) not in m.links:
        return pyo.Constraint.Skip
    return m.v_dot_sw_feed[p] == m.v_dot_link['SW', p]
model.sw_feed_link_coupling = pyo.Constraint(model.processes, rule=sw_feed_link_coupling_rule)


# solves for concentration of feed into ED
#BL - when I remove this constraint I get weird behavior. Dont know why but leaving it in
def ed_feed_concentration_rule(m, i):
    return m.f_ion_in['ED', i] == m.v_dot_in['ED'] * m.conc_feed['ED', i]
model.ed_feed_concentration = pyo.Constraint(model.ions, rule=ed_feed_concentration_rule)



 # lithium precipitation criteria is that concentrated ed stream needs 5 times as much sodium as lithium
# this has been changed from "equal" to "less than or equal to"
def ed_precipitation_ratio_rule(m):
    return m.concentration_dil['ED','Na'] == m.Ratio_ed_na_li * m.concentration_dil['ED','Li']
model.ed_precipitation_ratio = pyo.Constraint(rule=ed_precipitation_ratio_rule)


# lithium recovered is yield * what is entering ED
def li_recovery_rule(m):
    return m.li_recovered_ed == m.yield_ed_li * m.f_ion_in['ED','Li']
model.li_recovery = pyo.Constraint(rule=li_recovery_rule)


# connects lithium recovered to diluted stream
def ed_li_conc_mass_balance_rule(m):
    return m.li_recovered_ed == m.concentration_dil['ED','Li'] * m.v_dot_dil['ED']
model.ed_li_conc_mass_balance = pyo.Constraint(rule=ed_li_conc_mass_balance_rule)


# connects lithium recovered to what is still in diluted stream
def ed_li_dil_mass_balance_rule(m):
    return m.concentration_conc['ED','Li'] * m.v_dot_conc['ED'] == m.f_ion_in['ED','Li'] - m.li_recovered_ed
model.ed_li_dil_mass_balance = pyo.Constraint(rule=ed_li_dil_mass_balance_rule)



# this makes the concentration in both sides of the ED stream equal...not sure if this is the best choice 
# (although the concentrations would be equal entering the stream,
# but without setting some type of relationship it will make one side extremely concentrated)
#BL thinks that while this constraint does nto affect the results, it should be left in to keep Na stream tidy.
def ed_na_conc_mass_balance_rule(m):
    return m.concentration_conc['ED','Na'] == m.concentration_dil['ED', 'Na']
model.ed_na_conc_mass_balance = pyo.Constraint(rule=ed_na_conc_mass_balance_rule)


# end of new constraints


# amount of heat generated is a function of the enthalpy of the steam at that point and the extraction flow rate
def heat_generated_rule(m, c):
    return m.q_generated[c] == m.H_extract[c] * m.m_dot_extract_unit[c]
model.heat_generated = pyo.Constraint(model.cycles, rule=heat_generated_rule)

 
# amount of heat allocated to all processes from each power cycle configuration can not exceed the amount of heat generated from the power cycle
def heat_allocated_capacity_rule(m, c):
    return sum(m.q_allocated[c, p] for p in m.processes) == m.q_generated[c]
model.heat_allocated_capacity = pyo.Constraint(model.cycles, rule=heat_allocated_capacity_rule)


# amount of heat used by each process is the sum of the heat allocated to it from each cycle
def heat_used_balance_rule(m, p):
    return m.q_used[p] == sum(m.q_allocated[c, p] for c in m.cycles)
model.heat_used_balance = pyo.Constraint(model.processes, rule=heat_used_balance_rule)            


# heat cannot be allocated to a process unless the binary indicator for supplying heat from cycle c to process p is turned on
def heat_gating_rule(m, c, p):
    return m.q_allocated[c, p] <= m.Big_m_heat * m.y_q_source[c, p]
model.heat_gating = pyo.Constraint(model.cycles, model.processes, rule=heat_gating_rule)


# each process can only be supplied heat from one cycle configuration
#def single_heat_source_rule(m, p):
#    return sum(m.y_q_source[c, p] for c in m.cycles) <= 1
#model.single_heat_source = pyo.Constraint(model.processes, rule=single_heat_source_rule)


# general rule for heat used by each process, there are explicitly defined equations for ed and crystallization
def heat_used_general_rule(m, p):
    #if p == 'CRY':
    #    return pyo.Constraint.Skip
    return m.q_used[p] == m.Q_required[p] * m.v_dot_in[p]
model.heat_used_general = pyo.Constraint(model.processes, rule=heat_used_general_rule)

# electricity-related constraints

## electricity generated must either be sold or allocated to a process
def electricity_balance_rule(m, c):
    return m.elec_sold[c] == m.elec_generated_unit[c] - sum(m.elec_allocated[c, p] for p in m.processes)
model.electricity_balance = pyo.Constraint(model.cycles, rule=electricity_balance_rule)


# electricity used by each process is related to the electrical requirements of each process
def electricity_requirement_rule(m, p):
    return m.elec_used[p] == m.Elec_required[p] * m.v_dot_in[p]
model.electricity_requirement = pyo.Constraint(model.processes, rule=electricity_requirement_rule)


# electricity used by each process is the sum of the electricity allocated to it from each cycle
def electricity_used_balance_rule(m, p):
    return m.elec_used[p] == sum(m.elec_allocated[c, p] for c in m.cycles)
model.electricity_used_balance = pyo.Constraint(model.processes, rule=electricity_used_balance_rule)


#BL: this rule was removed because alocatede+sold=generated via electricity balance, so it is redundant
# amount of electricity allocated to all processes from each power cycle configuration can not exceed the amount of electricity generated from the power cycle
#def electricity_allocated_capacity_rule(m, c):
#    return sum(m.elec_allocated[c, p] for p in m.processes) <= m.elec_generated_unit[c]
#model.electricity_allocated_capacity = pyo.Constraint(model.cycles, rule=electricity_allocated_capacity_rule)



# mass conservation constraints

# seawater intake flow is limited by parameter, sums over all connections that include sw as a source node, however a single upstream and downstream node is also enforced in the model
def feed_intake_capacity_rule(m):
    return sum(m.v_dot_link['SW', p] for (upstream, p) in m.links if upstream == 'SW') <= m.V_dot_max
model.feed_intake_capacity = pyo.Constraint(rule=feed_intake_capacity_rule)


## mass flow in power cycle either has to be expanded in the turbines or extracted to provide heat to thermal processes, per cycle flow is limited by nominal conditions
def split_flow_rule(m, c):
    return m.m_dot_cycle_unit[c] + m.m_dot_extract_unit[c] == m.M_dot_max * model.y_cycle_active[c]
model.split_flow = pyo.Constraint(model.cycles, rule=split_flow_rule)


# conservation of volumetric flow through the inlet and outlet of each process
def process_flow_balance_rule(m, p):
    return m.v_dot_in[p] == m.v_dot_conc[p] + m.v_dot_dil[p]
model.process_flow_balance = pyo.Constraint(model.processes, rule=process_flow_balance_rule)


# ion conservation from inlet to concentrated and diluted stream
def ion_conservation_rule(m, p, i):
    if p == 'CRY':
        return pyo.Constraint.Skip
    else:
        return m.f_ion_in[p, i] == m.v_dot_conc[p] * m.concentration_conc[p, i] + m.v_dot_dil[p] * m.concentration_dil[p, i]
model.ion_conservation = pyo.Constraint(model.processes, model.ions, rule=ion_conservation_rule)


# the volumetric flow into a process is the sum of the flow on the links connected to the process upstream
def inlet_flow_balance_rule(m, p):
    inflow_from_links = sum(m.v_dot_link[q, downstream] for (q, downstream) in m.links if downstream == p)
    return m.v_dot_in[p] == inflow_from_links
model.inlet_flow_balance = pyo.Constraint(model.processes, rule=inlet_flow_balance_rule)



# logical constraints

# processes must be turned on to have a concentration
# def logical_concentration_inactive_rule(m, p, i):
#     return m.concentration_conc[p,i] <= m.Little_m + m.Conc_max[p, i] * m.y_process_active[p]
# model.logical_concentration_inactive = pyo.Constraint(model.processes, model.ions, rule=logical_concentration_inactive_rule)


# connects the electricity generated to the piecewise variable elec_generated_unmasked through the binary variable activation
def logical_mask_cycle_generation_rule(m, c):
    return m.elec_generated_unit[c] == m.elec_generated_unmasked[c] * m.y_cycle_active[c]
model.logical_mask_cycle_generation = pyo.Constraint(model.cycles, rule=logical_mask_cycle_generation_rule)


# if a link is turned on there has to be flow on it
def logical_no_flow_means_inactive_rule(m, q, p):
    return m.v_dot_link[q, p] >= m.Little_m * m.y_link_active[q, p]
model.logical_no_flow_means_inactive = pyo.Constraint(model.links, rule=logical_no_flow_means_inactive_rule)


# process can only have inlet flow if its binary activation is turned on
def logical_inflow_maximum_activation_rule(m, p):
    return m.v_dot_in[p] <= m.Big_m * m.y_process_active[p]
model.logical_inflow_maximum_activation = pyo.Constraint(model.processes, rule=logical_inflow_maximum_activation_rule)


# process must have inlet flow if it's binary activation is turned on
def logical_inflow_minimum_activation_rule(m, p):
    return m.v_dot_in[p] >= m.Little_m * m.y_process_active[p]
model.logical_inflow_minimum_activation = pyo.Constraint(model.processes, rule=logical_inflow_minimum_activation_rule)
   

# seawater feed can only go to one process to start the chain of processes
def logical_single_feed_target_rule(m):
    return sum(m.y_link_active['SW', p] for p in m.processes if ('SW', p) in m.links) <= 1
model.logical_single_feed_target = pyo.Constraint(rule=logical_single_feed_target_rule)

# each downstream process can be fed by only one upstream process
def logical_single_upstream_rule(m, p):
    return sum(m.y_link_active[upstream, p] for (upstream, downstream) in m.links if downstream==p) <= 1
model.logical_single_upstream = pyo.Constraint(model.processes, rule=logical_single_upstream_rule)


# each process can send outlet flow to only one downstream process, checks to make sure there is a downstream process and if not skips the constraint
def logical_single_downstream_rule(m, p):
    terms = [m.y_link_active[p, downstream] for (upstream, downstream) in m.links if upstream == p]
    if not terms:
        return pyo.Constraint.Skip
    return sum(terms) <= 1
model.logical_single_downstream = pyo.Constraint(model.processes, rule=logical_single_downstream_rule)


# link between processes can only be activated if the upstream process is activated
def logical_link_activation_upstream_rule(m, upstream, downstream):
    if upstream == 'SW':
        return pyo.Constraint.Skip
    return m.y_link_active[upstream, downstream] <= m.y_process_active[upstream]
model.logical_link_activation_upstream = pyo.Constraint(model.links, rule=logical_link_activation_upstream_rule)


# link between processes can only be activated if the downstream process is activated
def logical_link_activation_downstream_rule(m, upstream, downstream):
    return m.y_link_active[upstream, downstream] <= m.y_process_active[downstream]
model.logical_link_activation_downstream = pyo.Constraint(model.links, rule=logical_link_activation_downstream_rule)


# the flow from one link to another is 0 if the link is inactive; otherwise it is limited by a big M
def logical_link_capacity_rule(m, upstream, downstream):
    return m.v_dot_link[upstream, downstream] <= m.Big_m * m.y_link_active[upstream, downstream]
model.logical_link_capacity = pyo.Constraint(model.links, rule=logical_link_capacity_rule)



# the ion flow to each process is related to what is feeding it 
def logical_ion_mass_in_rule(m, p, i):
    
    seawater_term = 0.0
    if ('SW', p) in m.links:
        seawater_term = m.v_dot_sw_feed[p] * m.Conc_sw[i] * m.y_link_active['SW', p]
           
    upstream_term = sum(
        m.v_dot_link[upstream, p] * m.concentration_conc[upstream, i]
        for (upstream, downstream) in m.links if downstream == p and upstream != 'SW'
    )
    
    return m.f_ion_in[p, i] == seawater_term + upstream_term
model.logical_ion_mass_in = pyo.Constraint(model.processes, model.ions, rule=logical_ion_mass_in_rule)

""" BL - after rerunning tests, I concur this is redundant and have removed it again
# flow routing depending on upstream process
def logical_inlet_flow_routing_rule(m, p):
    
    seawater_term = 0.0
    
    if ('SW', p) in m.links:
        seawater_term = m.v_dot_sw_feed[p] * m.y_link_active['SW', p]
    
    upstream_term = sum(
            m.v_dot_conc[upstream] * m.y_link_active[upstream, p]
            for (upstream, downstream) in m.links if downstream == p and upstream != 'SW')
    
    return m.v_dot_in[p] == seawater_term + upstream_term

model.logical_inlet_flow_routing = pyo.Constraint(model.processes, rule=logical_inlet_flow_routing_rule)
"""

# ensures that the flow on the links is correct depending on the processes being used
def logical_link_flow_balance(m, upstream, p):
    
    if (upstream, p) not in m.links or upstream == 'SW':
        return pyo.Constraint.Skip

    else:
        return m.v_dot_link[upstream, p] == m.v_dot_conc[upstream] * m.y_link_active[upstream, p]
    
model.logical_link_flow_balance = pyo.Constraint(model.links, rule=logical_link_flow_balance)

# ro specific constraints

# diluted (potable) water is a function of recovery ratio
def recovery_ratio_ro_rule(m):
    return m.v_dot_dil['RO'] == m.Recovery_ro * m.v_dot_in['RO']
model.recovery_ratio_ro = pyo.Constraint(rule=recovery_ratio_ro_rule)


# no ions in diluted ro stream
def diluate_concentration_ro_rule(m, p, i):
    if p!= 'RO':
        return pyo.Constraint.Skip
    return m.concentration_dil['RO', i] == 0.0
model.diluate_concentration_ro = pyo.Constraint(model.processes, model.ions, rule=diluate_concentration_ro_rule)



# med specific constraints

# diluted water is a function of recovery ratio - recovery ratio for med depends on if it is fed by seawater or ro brine
def logical_recovery_ratio_med_rule(m):
    return m.v_dot_dil['MED'] == (
        m.v_dot_in['MED'] * m.Recovery_sw_to_med * (1 - m.y_link_active['RO','MED']) + 
        m.v_dot_in['MED'] * m.Recovery_ro_to_med * m.y_link_active['RO','MED'])
    
model.logical_recovery_ratio_med = pyo.Constraint(rule=logical_recovery_ratio_med_rule)


# no ions in diluted med stream
def diluate_concentration_med_rule(m, p, i):
    if p!= 'MED':
        return pyo.Constraint.Skip
    return m.concentration_dil['MED', i] == 0.0
model.diluate_concentration_med = pyo.Constraint(model.processes, model.ions, rule=diluate_concentration_med_rule)



# ed specific constraints

# the mass concentration of lithium in the concentrated ed stream is limited by optimal concentration for precipitation
# def ed_li_maximum_concentration_rule(m):
#     return m.concentration_conc['ED', 'Li'] <= m.Conc_max_ed_li
# model.ed_li_maximum_concentration = pyo.Constraint(rule=ed_li_maximum_concentration_rule)


# the mass concentration of sodium in the concentrated ed stream is limited by optimal concentration for precipitation
# def ed_na_maximum_concentration_rule(m):
#     return m.concentration_conc['ED', 'Na'] <= m.Conc_max_ed_na
# model.ed_na_maximum_concentration = pyo.Constraint(rule=ed_na_maximum_concentration_rule)


# chloride ions are split between the two streams
def ed_cl_concentration_rule(m):
    return m.concentration_conc['ED', 'Cl'] == m.concentration_dil['ED', 'Cl']
model.ed_cl_concentration = pyo.Constraint(rule=ed_cl_concentration_rule)



# crystallization constraints

# there is no concentrated brine from crystallization, only potable water and solid salts
def cry_no_concentrated_flow_rule(m):
    return m.v_dot_conc['CRY'] == 0.0
model.cry_no_concentrated_flow = pyo.Constraint(rule=cry_no_concentrated_flow_rule)

    
# there are no ions in the diluted stream outflowing from crystallization
def diluate_concentration_cry_rule(m, i):
    return m.concentration_dil['CRY', i] == 0.0
model.diluate_concentration_cry = pyo.Constraint(model.ions, rule=diluate_concentration_cry_rule)

# OBJECTIVE FUNCTION

profit = (
    model.Price_li * 8760 * model.li_recovered_ed + 
    model.Price_water * 8760 * (model.v_dot_dil['RO'] + model.v_dot_dil['MED'] + model.v_dot_dil['CRY']) +
    model.Price_elec * 8760 * sum(model.elec_sold[c] for c in model.cycles) +
    model.Price_salt * 8760 * sum(model.f_ion_in['CRY', i] for i in model.ions))

cost   = (sum(model.K_cycle * model.elec_generated_unit[c] for c in model.cycles)/30 +
         sum(model.K_process[p] * model.v_dot_in[p] * 24 for p in model.processes)/30
)
model.obj = pyo.Objective(expr=(profit-cost), sense=pyo.maximize)


solver = pyo.SolverFactory("gurobi")
solver.options['NonConvex'] = 2
results = solver.solve(model, tee=False)

# RESULTS 
def safe_val(v, fmt="{:.2f}"):
    """Return variable value or 'not solved' if uninitialized."""
    if v.value is None:
        return "not solved"
    try:
        return fmt.format(pyo.value(v))
    except:
        return v.value

def nval(x):
    try:
        v = pyo.value(x)
        return None if v is None else float(v)
    except Exception:
        return None

def fmt_flow(v):
    return "not solved" if v is None else f"{v:.2f} m³/h"

def fmt_mass(v):
    return "not solved" if v is None else f"{v:.2f} kg/h"

def fmt_conc(v):
    if v is None:
        return "not solved"
    mgL = v * 1000.0
    return f"{v:.6f} kg/m³ ({mgL:.1f} mg/L)"

def fmt_resid(v):
    if v is None:
        return "n/a"
    tag = "" if abs(v) <= 1e-6 else "  <-- check"
    return f"{v:+.6e} kg/h{tag}"


print(f" li_yield = {safe_val(model.yield_ed_li)}")


print("\n================== POWER CYCLE RESULTS ==================")
for c in model.cycles:
    print(f"\nCycle {c}:")
    print(f"  ṁ_cycle_unit        = {safe_val(model.m_dot_cycle_unit[c])} kg/s")
    print(f"  ṁ_extract_unit      = {safe_val(model.m_dot_extract_unit[c])} kg/s")
    print(f"  Ẇ_generated_unit    = {safe_val(model.elec_generated_unit[c])} kWe")
    print(f"  Ẇ_sold_unit         = {safe_val(model.elec_sold[c])} kWe")
    print(f"  Q_generated         = {safe_val(model.q_generated[c])} kWth")


print("\n================== PROCESS RESULTS ==================")
for p in model.processes:
    print(f"\nProcess {p}:")
    print(f"  Active?             = {safe_val(model.y_process_active[p], '{:.0f}')}")
    print(f"  v̇_in                = {safe_val(model.v_dot_in[p])} m³/h")
    print(f"  v̇_conc              = {safe_val(model.v_dot_conc[p])} m³/h")
    print(f"  v̇_dil               = {safe_val(model.v_dot_dil[p])} m³/h")
    print(f"  Electricity used    = {safe_val(model.elec_used[p])} kWe")
    print(f"  Heat used           = {safe_val(model.q_used[p])} kW-th")

    for i in model.ions:
        print(f"    Ion {i}:")
        print(f"      f_in            = {safe_val(model.f_ion_in[p, i])} kg/h")
        print(f"      c_conc          = {fmt_conc(nval(model.concentration_conc[p, i]))}")
        print(f"      c_dil           = {fmt_conc(nval(model.concentration_dil[p, i]))}")


print("\n================== ACTIVE LINKS ==================")
active = False
for (q, p) in model.links:
    y = nval(model.y_link_active[q, p])
    v = nval(model.v_dot_link[q, p])
    if y is not None and (y > 0.5 or (v is not None and v > 1e-3)):
        print(f"  {q} → {p}:  y_active={int(round(y))}, v̇ = {v:.2f} m³/h")
        active = True
if not active:
    print("  (No active links)")


print("\n================== OBJECTIVE VALUE ==================")
try:
    print(f"  Profit = {pyo.value(model.obj):,.2f} $")
except:
    print("  Profit = not solved")
