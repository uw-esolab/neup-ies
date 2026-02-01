import pyomo.environ as pyomo
import pandas as pd


from power_cycle_calcs import compute_power_cycle_nominals


# scalar parameters

def add_scalar_params(model, sweep_factor=1.0):
    
    


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
model.M_lt_tes_init           = pyomo.Param(initialize=0)                                                       # units: kg,           initial storage inventory of low-temperature tes
model.M_ht_tes_init           = pyomo.Param(initialize=0)                                                       # units: kg,           initial storage inventory of high-temperature tes 
model.Cost_ramp_power         = pyomo.Param(initialize=43.75*W_dot_gen/M_dot_steam_cyc_nom)                     # units: $/ 
model.Cost_ht_tes             = pyomo.Param(initialize=Price_storage_salt/model.Convert_time.value*Cp_salt*Delta_T_salt)     # units: $/kg      cost of high-temperature (molten salt) storage
model.Cost_lt_tes             = pyomo.Param(initialize=Price_storage_water/model.Convert_time.value*Delta_H_water)           # units: $/kg      cost of low-temperature (pressurized water) storage
model.Price_water             = pyomo.Param(initialize=0.0021)                                                  # units: $/kg         selling price of water 
model.Price_elec              = pyomo.Param(model.T, initialize={t: float(Schedule_elec[t-1]/1000) for t in model.T})            # units: $/kWh       selling price of electricity(divide by 1000 to get from MWh to kWh)
model.Slope_power_curve       = pyomo.Param(initialize=1.075)
model.Interc_power_curve      = pyomo.Param(initialize=0.0748)
model.F_lpt_post              = pyomo.Param(initialize=0.186)                                                   # units: --,          fraction of total LPT power generated downstream of extraction point



model.Recovery_ro             = pyomo.Param(initialize=0.50)   # units: --,     recovery ratio for ro fed seawater
model.Recovery_sw_to_med      = pyomo.Param(initialize=0.77)   # units: --,     recovery ratio for med fed seawater
model.Recovery_ro_to_med      = pyomo.Param(initialize=0.53)   # units: --,     recovery ratio for med fed ro concentrate
model.Price_li                = pyomo.Param(initialize=70.0) # units: $/kg,   sales price of lithium
model.V_dot_min                = pyomo.Param(initialize=1e-04)  # units: --,     little m used for logical constraints
model.Ratio_ed_na_li          = pyomo.Param(initialize=5.0)    # units: --,     optimal ratio between sodium and lithium in concentrated stream leaving ed for precipitation            
model.V_dot_max               = pyomo.Param(initialize=40000)  # units; m3/h,   maximum seawater inflow for desalination, based on current largest desalination plant in world 
model.Conc_sw                 = pyomo.Param(model.ions, initialize={'Li': 0.00018, 'Na': 10.8, 'Cl': 19.3})                                                               # units: kg/m3     mass concentration of ions in seawater
model.Elec_required           = pyomo.Param(model.processes, initialize={'RO': 4.5, 'MED': 1.0, 'ED': 3.5,'CRY': 3.5})                                                # units: kWe/m3    electrical energy required per unit of feed 
model.Q_required              = pyomo.Param(model.processes, initialize={'RO': 0.0,  'MED': 25.0 ,'ED': 0.0, 'CRY': 90})                                                 # units: kW-th/m3  heat required per unit of feed (for ed and crystallization these values are calculated based on equations in the model)
model.K_process               = pyomo.Param(model.processes, initialize={'RO': 500.0,'MED': 900.0,'ED': 200.0, 'CRY': 950.0})                                                 # units: $/m3/day  capex for each process
model.yield_ed_li             = pyomo.Param(initialize=1.0)    # units: --,   li recovery fraction, assumed 100%

