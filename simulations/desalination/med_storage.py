import pyomo.environ as pyomo
import matplotlib.pyplot as plt
import random
import csv

time_steps = 48
Vmax = 7500                                     #m^3/hr
percent_reduction_power = 0.2
maximum_discharge = Vmax
V0 = 0.5 * Vmax
Ymax_MW = 500   

def create_electric_price_schedule(time_steps):
    electric_price_schedule = {}
    for i in range(time_steps):
        electric_price_schedule[i] = random.uniform(40, 100)
    return electric_price_schedule

def create_water_price_schedule(time_steps):
    water_price_schedule = {}
    for i in range(time_steps):
        water_price_schedule[i] = random.uniform(1, 3)
    return water_price_schedule

def create_inflow_schedule(time_steps):
    inflow_schedule = {}
    for i in range(time_steps):
        inflow_schedule[i] = 0.5 * Vmax
    return inflow_schedule


model = pyomo.ConcreteModel()

# parameters
 
# number of time steps
model.nt = pyomo.Param(initialize = time_steps, domain = pyomo.Integers)
# set of time steps
model.T = pyomo.Set(initialize = range(model.nt()))
# sales price of discharge at each time step
model.price_water = pyomo.Param(model.T, initialize = create_water_price_schedule(time_steps))
# sales price of electricty at each time step
model.price_electric = pyomo.Param(model.T, initialize = create_electric_price_schedule(time_steps))
# water added at each time step
model.inflow = pyomo.Param(model.T, initialize = create_inflow_schedule(time_steps))
 # maximum size of storage tank
model.Vmax = pyomo.Param(initialize = Vmax)
# initial volume of water in tank
model.V0 = pyomo.Param(initialize = V0)
# maximum discharge from tank
model.dmax = pyomo.Param(initialize = maximum_discharge)
# maximum output from power-producing system
model.Ymax = pyomo.Param(initialize = Ymax_MW)



# variables

# volume of water stored in tank - dependent variable
model.V = pyomo.Var(model.T, domain = pyomo.NonNegativeReals)
# water discharged from tank
model.d = pyomo.Var(model.T, domain = pyomo.NonNegativeReals)
# water charged into storage tank
model.c = pyomo.Var(model.T, domain = pyomo.NonNegativeReals)
# power produced by turbine - dependent variable
model.y = pyomo.Var(model.T, domain = pyomo.NonNegativeReals)
# amount of water sent to power-producing system
#model.r = pyomo.Var(model.T, domain = pyomo.NonNegativeReals)

# objective function

# sum the discharge times price for all time steps 
def objective_func(model):
    return sum((model.y[t] * model.price_electric[t]) + (model.d[t] * model.price_water[t]) for t in model.T)
model.objective = pyomo.Objective(rule = objective_func, sense = pyomo.maximize)


# constraints

# conservation of mass on system - water going to power-producing side & MED storage tank
def constr_mass_balance_system(model, t):
    return model.inflow[t] >= model.c[t]
model.constr_mass_balance_system = pyomo.Constraint(model.T,rule = constr_mass_balance_system)

# electricy discharged limited by size of turbine
def constr_electric_capacity(model, t):
    return model.y[t] <= model.Ymax
model.constr_electric_capacity = pyomo.Constraint(model.T, rule = constr_electric_capacity)

# water discharged limited by size of MED
def constr_discharge_capacity(model, t):
    return model.d[t] <= model.dmax
model.constr_discharge_capacity = pyomo.Constraint(model.T, rule = constr_discharge_capacity)

# water stored is limited by size of storage tank
def constr_water_storage(model, t):
    return model.V[t] <= model.Vmax
model.constr_water_storage = pyomo.Constraint(model.T, rule = constr_water_storage)

# mass balance on storage tank
def constr_storage_tank_balance(model, t):
    if t == 0:
        return model.V[t] == model.V0 - model.d[t] + model.c[t]
    else:
        return model.V[t] == model.V[t-1] - model.d[t] + model.c[t]
model.constr_storage_tank_balance = pyomo.Constraint(model.T, rule = constr_storage_tank_balance)

# power production as a function of steam flow to turbine
# dependent on the fraction of the steam sent to the turbine which is 1 - fraction of steam charged to MED storage
def constr_power_production(model, t):
    return model.y[t] <= (1 - percent_reduction_power * (model.c[t]/model.inflow[t])) * model.Ymax
model.constr_power_production = pyomo.Constraint(model.T, rule = constr_power_production)

solver = pyomo.SolverFactory('gurobi')

results = solver.solve(model, keepfiles = False, logfile = 'solve.log')

print(model.display())
print(f"time\tprice_water\tMED_charge\tMED_discharge\ttank_storage\tprice_elec\tpower_sold")
for t in model.T:
    print(f'{t}\t{model.price_water[t]:.2f}\t{model.c[t]():>5.1f}\t{model.d[t]():>5.1f}\t{model.V[t]():>5.1f}\t{model.price_electric[t]:.2f}\t{model.y[t]():>5.1f}')

