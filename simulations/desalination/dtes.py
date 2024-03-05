"This model generates data for different low-temperature TES prices ranging from 10% base case to 200% base case."
import pyomo.environ as pyomo
import pandas as pd
import numpy as np

# ----------------------------
# prep

counter_dtes                                      = 1

W_dot_nom                                         = 494 
cost_dtes                                         = 20
sens_dtes                                         = [0.1*cost_dtes, 0.2*cost_dtes, 0.3*cost_dtes, 0.4*cost_dtes, 0.5*cost_dtes, 0.6*cost_dtes, 0.7*cost_dtes, 0.8*cost_dtes, 0.9*cost_dtes, 1.0*cost_dtes, 1.1*cost_dtes, 1.2*cost_dtes, 1.3*cost_dtes, 1.4*cost_dtes, 1.5*cost_dtes, 1.6*cost_dtes, 1.7*cost_dtes, 1.8*cost_dtes, 1.9*cost_dtes, 2.0*cost_dtes]



for dtes in sens_dtes:
    
        
    N                                             = 8760
    
    # Some preliminary sizing calculations
    Q_dot_LFR                                     = 950                                                                               #[MW]
    eta_cycle_ref                                 = 0.52
    W_dot_gen                                     = W_dot_nom
    cp_salt                                       = 1.5                                                                               #[kJ/kg-K]
    DeltaT_ctes                                   = 565-330                                                                           #[K]
    DeltaH_csteam                                 = 2126                                                                              #[kJ/kg]           (h(150bar,555C) - h(150bar,300C))
    DeltaH_dsteam                                 = 2494                                                                              #[kJ/kg]           (h(4.8bar,x=.98) - h(4.5bar,50C))  # steam extracted from turbine
    DeltaH_dtes                                   = 421.3                                                                             #[kJ/kg]           (h(5bar,140C) - h(5bar,40C))       # water through HX to desal hot storage 
    
    
    # price of desalinated water ($/kg)
    price_dist                                    = 0.0021 
                                                                             
    # cost of high-temperature storage ($/kWh)
    cost_ctes                                     = 29.8  

    
    m_dot_salt                                    = Q_dot_LFR*1e3 /(cp_salt * DeltaT_ctes)                                           #[kg/s]
    m_dot_csteam                                  = W_dot_gen/eta_cycle_ref*1e6 /(DeltaH_csteam*1e3)                                 #[kg/s]
    
    f_lpt                                         = 0.53                                                                                                 #nominal fraction of steam flow that makes it to LPT stage 3
    
    m_dot_esteam                                  = m_dot_csteam*f_lpt                                                                                   #nominal extraction flow rate
    m_dot_dtes                                    = m_dot_esteam*DeltaH_dsteam/DeltaH_dtes                                                               #nominal flow rate of water into the desal storage tank
    
    
    sched_inflow                                  = [m_dot_salt for n in range(N)]                                                                       #LFR salt into storage
   
    sched_eff      = pd.read_csv('../../../input/sched_eff.csv', usecols=['Efficiency'])
    sched_eff      = sched_eff.values

    sched_elec     = pd.read_csv('../../../input/DIABLOCN_2_N001_PRICES.csv', usecols=['LMP'])
    sched_elec     = sched_elec.values
    
    # # -------------------------------
    model = pyomo.ConcreteModel()
    
    # #  ==================================================
    # #       parameters
    # #  ==================================================
    model.nt                                      = pyomo.Param(initialize = len(sched_inflow), domain = pyomo.Integers)                 #[-]             Number of time steps
    model.T                                       = pyomo.Set(initialize = range(1,N+1))                                                #[-]              Set of time steps
    model.Delta_t                                 = pyomo.Param(initialize = 1)                                                         #[hr]             Time step duration
    model.M_dot_LFR                               = pyomo.Param(model.T, initialize = dict(zip(model.T, sched_inflow)))                 #[kg/s]           Mass flow produced by the LFR at time t
    model.eta                                     = pyomo.Param(model.T, initialize = dict(zip(model.T, sched_eff)))                    #[-]              Ambient temperature efficiency modifier in the power cycle at time t
    model.P_elec                                  = pyomo.Param(model.T, initialize = dict(zip(model.T, sched_elec)))                   #[$/MWh]          Price at which electricity is sold at time t
    model.P_dist                                  = pyomo.Param(initialize=price_dist)                                                  #[$/kg]           Price at which distillate is sold at time t
    model.C_c_ramp                                = pyomo.Param(initialize=43.75*W_dot_gen/m_dot_csteam)                                #[$/Delta kg/s]   Cost associated with change in electric power production, from Gabriel's paper
    model.C_d_ramp                                = pyomo.Param(initialize=0.2)                                                         #[$/Delta kg/s]   Cost associated with change in distillate production 
    model.K_chx                                   = pyomo.Param(initialize=m_dot_csteam/m_dot_salt)                                     #[(kg/s)/(kg/s)]  Conversion constant for salt mass flow to steam mass flow in the cycle heat exchanger
    model.K_chpt                                  = pyomo.Param(initialize=W_dot_gen/m_dot_csteam)                                      #[MW/(kg/s)]      Conversion constant for HPT steam mass flow to power in the cycle
    model.K_clpt                                  = pyomo.Param(initialize=W_dot_gen/m_dot_esteam)                                      #[MW/(kg/s)]      Conversion constant for LPT extraction steam mass flow to lost power in the cycle
    model.K_MED                                   = pyomo.Param(initialize=1095.3)                                                      #[$/m^3/day]      Slope of MED capital cost graph
    model.K_dint                                  = pyomo.Param(initialize=2.05E07)                                                     #[$]              Constant equation intercept for desalination capital costs
    model.alpha                                   = pyomo.Param(initialize=1.075)                                                       #[]               Slope of power curve
    model.beta                                    = pyomo.Param(initialize=0.186)                                                       #[]               Fraction of power generated in the low-pressure turbine past the extraction point
    model.K_cint                                  = pyomo.Param(initialize=0.0748)                                                      #[-]              Constant equation intercept for cycle mass flow to power conversion
    model.K_vdot                                  = pyomo.Param(initialize=86.4)                                                        #[m^3/day/kg/s]   Conversion between m^3/day and kg/s for distillate production
    model.K_t                                     = pyomo.Param(initialize=3600)                                                        #[s/h]            Converstion constant for time
    model.F_lpt                                   = pyomo.Param(initialize=0.53)                                                        #[-]              Max fraction of steam mass flow available for extraction
    model.W_dot_ramp_max                          = pyomo.Param(initialize=0.1*W_dot_gen)                                               #[Delta MW/hr]    Maximum allowable power system ramp rate 
    model.K_dhx                                   = pyomo.Param(initialize=m_dot_dtes/m_dot_esteam)                                     #[(kg/s)/(kg/s)]  Conversion constant for HPT extraction steam mass flow to desal storage charge mass flow
    model.M_cm_min                                = pyomo.Param(initialize=0)                                                           #[kg]             Cycle thermal storage minimum inventory
    model.M_dm_min                                = pyomo.Param(initialize=0)                                                           #[kg]             Desal thermal storage minimum inventory
    model.K_d                                     = pyomo.Param(initialize=0.70)                                                        #[(kg/s)/(kg/s)]  Rate of distillate production per unit mass flow into the desal system
    model.M_dot_ch_init                           = pyomo.Param(initialize=0.0)                                                         #[kg]             Initial cycle thermal storage inventory
    model.M_dot_dh_init                           = pyomo.Param(initialize=0.0)                                                         #[kg]             Initial desal thermal storage inventory
    model.W_dot_max                               = pyomo.Param(initialize=W_dot_gen)                                                   #[MW]             Maximum power produced at each time step
    model.W_dot_min                               = pyomo.Param(initialize=W_dot_gen*0.25)                                              #[MW]             Minimum power produced at each time step
    model.V_dot_min                               = pyomo.Param(initialize=0)                                                           #[kg/s]           Minimum amount of distillate produced at each time step
    model.C_cs                                    = pyomo.Param(initialize=cost_ctes)                                                   #[$/kWh]          Cost associated with cycle thermal energy storage
    model.C_ds                                    = pyomo.Param(initialize=dtes)                                                        #[$/kWh]          Cost associated with desal thermal energy storage
    model.T_a                                     = pyomo.Param(initialize=30)                                                          #[yr]             Amortization period for desal and thermal storage subsystems
    
    
    #  ==================================================
    #       Variables
    #  ==================================================
    model.m_ch                                    = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg]             Inventory in cycle-side hot storage at time step t
    model.m_dot_cs                                = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg/s]           Mass flow drawn from cycle-side hot storage at time step t
    model.m_dot_hpt                               = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg/s]           Mass flow produced in the cycle at the high pressure turbine inlet at time step t
    model.m_dot_hpt_max                           = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg/s]           Maximum mass flow to the high pressure turbine
    model.m_dot_e                                 = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg/s]           Mass flow produced extracted and sent to the desal heat exchanger at time step t
    model.m_dh                                    = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg]             Inventory in the desal-side hot storage at time step t
    model.m_dot_ds                                = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg/s]           Mass flow consumed by the desal system at time step t
    model.w_dot                                   = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[MW]             Electric power produced by the turbine system at time step 
    model.v_dot                                   = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                   #[kg/s]           Distillate produced by the desal system at time step t
    model.w_dot_Delta_up                          = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)                               #[MW]             Positive change in electric power produced relative to the previous time step
    model.w_dot_Delta_dn                          = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)                               #[MW]             Negative change in electric power produced relative to the previous time step
    model.v_dot_Delta_up                          = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)                               #[kg/s]           Positive change in distillate production relative to the previous time step
    model.v_dot_Delta_dn                          = pyomo.Var(model.T-[1], domain=pyomo.NonNegativeReals)                               #[kg/s]           Negative change in distillate production relative to the previous time step
    model.M_cm_max                                = pyomo.Var(domain=pyomo.NonNegativeReals)                                            #[kWh]            Cycle thermal storage maximum inventory
    model.M_dm_max                                = pyomo.Var(domain=pyomo.NonNegativeReals)                                            #[kWh]            Desal thermal storage maximum inventory
    model.V_dot_max                               = pyomo.Var(domain=pyomo.NonNegativeReals)                                            #[kg/s]           Size of desalination system
    model.C_desal                                 = pyomo.Var(domain=pyomo.NonNegativeReals)                                            #[$]              Cost of desalination system
    model.hpt_Delta_up                            = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                            #[kg/s]
    model.hpt_Delta_dn                            = pyomo.Var(model.T, domain=pyomo.NonNegativeReals)                                            #[kg/s]
    
    
    #  ==================================================
    #   Objective
    def objective(model):
        # return sum(model.w_dot[t] for t in model.T)
        return sum([model.P_elec[t]*model.w_dot[t]*model.Delta_t + model.P_dist*model.v_dot[t]*model.K_t*model.Delta_t for t in model.T]) \
            - sum([model.C_c_ramp*(model.hpt_Delta_up[t]+model.hpt_Delta_dn[t]) for t in (model.T-[1])]) \
            - sum([model.C_d_ramp*(model.v_dot_Delta_up[t]+model.v_dot_Delta_dn[t]) for t in (model.T-[1])]) \
                - (model.C_cs*model.M_cm_max)/model.T_a - (model.C_ds*model.M_dm_max)/model.T_a - model.C_desal/model.T_a
    model.objective = pyomo.Objective(rule=objective, sense=pyomo.maximize)
    
    #  ==================================================
    #   Constraints
    #  ==================================================
    
    def constr_C_desal(model, t):
        return model.C_desal                     == model.K_MED*(model.V_dot_max*model.K_vdot) + model.K_dint
    model.constr_C_desal                         =  pyomo.Constraint(model.T, rule=constr_C_desal)
    
    
    
    # ------- Power cycle subsystem ------------- 
    def constr_csmass(model,t):
        if t>1: 
            return model.m_ch[t]                 == (model.M_dot_LFR[t] - model.m_dot_cs[t])*model.K_t*model.Delta_t + model.m_ch[t-1] 
        else:
            return model.m_ch[t]                 == (model.M_dot_LFR[t] - model.m_dot_cs[t])*model.K_t*model.Delta_t + model.M_dot_ch_init
    model.constr_csmass                          =  pyomo.Constraint(model.T, rule=constr_csmass)
    
    #Inventory limits on the cycle storage 
    def constr_csmmax(model,t):
        return model.m_ch[t]                     <= (model.M_cm_max*model.K_t)/(cp_salt*DeltaT_ctes)
    model.constr_csmmax                          =  pyomo.Constraint(model.T, rule=constr_csmmax)
    def constr_csmmin(model,t):
        return model.m_ch[t]                     >= model.M_cm_min
    model.constr_csmmin                          =  pyomo.Constraint(model.T, rule=constr_csmmin)
    
    
    #Mass flow conversion from cycle storage to HPT inlet
    def constr_chx(model,t):
        return model.m_dot_hpt[t]                == model.K_chx * model.m_dot_cs[t]
    model.constr_chx                             =  pyomo.Constraint(model.T, rule=constr_chx)
    
    
    # Conversion of mass flow to power in the cycle
    def constr_cmass(model, t):
        return model.w_dot[t]                    == ((model.eta[t]*W_dot_gen)/m_dot_csteam) * (model.alpha*model.m_dot_hpt[t] - model.beta*model.m_dot_e[t] - model.K_cint*m_dot_csteam)
    model.constr_cmass                           =  pyomo.Constraint(model.T, rule=constr_cmass)
                                                             
    
    
    # Mass flow limit on the extraction
    def constr_mext(model, t):
        return model.m_dot_e[t]                  <= model.F_lpt * model.m_dot_hpt[t]
    model.constr_mext                            =  pyomo.Constraint(model.T, rule=constr_mext)
    
    
    # Range limits on power
    def constr_wbounds_up(model, t):
        return model.m_dot_hpt[t]                <= model.m_dot_hpt_max[t]
    model.constr_wbounds_up                      =  pyomo.Constraint(model.T, rule=constr_wbounds_up)
    
    
    def constr_wbounds_dn(model,t):
        return model.m_dot_hpt_max[t]            >= model.W_dot_min/model.K_chpt
    model.constr_wbounds_dn                      =  pyomo.Constraint(model.T, rule=constr_wbounds_dn)
    
    
    # Conversion between maximum power and mass flow to hpt
    def constr_hpt_flow(model,t):
        return model.m_dot_hpt_max[t]            == model.W_dot_max * model.eta[t]/model.K_chpt
    model.constr_hpt_flow                        =  pyomo.Constraint(model.T, rule=constr_hpt_flow)
    
    
    # Power ramping tracking
    def constr_hpt_track_up(model, t):
        return model.hpt_Delta_up[t]             >= model.m_dot_hpt[t] - model.m_dot_hpt[t-1]
    model.constr_hpt_track_up                    =  pyomo.Constraint(model.T-[1], rule=constr_hpt_track_up)
    
    def constr_hpt_track_dn(model, t):
        return model.hpt_Delta_dn[t]             >= model.m_dot_hpt[t-1] - model.m_dot_hpt[t]
    model.constr_hpt_track_dn                    =  pyomo.Constraint(model.T-[1], rule=constr_hpt_track_dn)
    
    
    # Power ramping limits
    def constr_cramp_up(model,t):
        return model.hpt_Delta_up[t]             <= model.W_dot_ramp_max/model.K_chpt * model.Delta_t
    model.constr_cramp_up                        =  pyomo.Constraint(model.T-[1], rule=constr_cramp_up)
    
    def constr_cramp_dn(model,t):
        return model.hpt_Delta_dn[t]             <= model.W_dot_ramp_max/model.K_chpt * model.Delta_t
    model.constr_cramp_dn                        =  pyomo.Constraint(model.T-[1], rule=constr_cramp_dn)  
    
    
    
    # Range limits on power
    # def constr_wbounds_up(model, t):
    #     return model.w_dot[t]                    <= model.W_dot_max * model.eta[t]
    # model.constr_wbounds_up                      =  pyomo.Constraint(model.T, rule=constr_wbounds_up)
    # def constr_wbounds_dn(model, t):
    #     return model.w_dot[t]                    >= model.W_dot_min
    # model.constr_wbounds_dn                      =  pyomo.Constraint(model.T, rule=constr_wbounds_dn)
    
    
    # # Power ramping tracking
    # def constr_wtrack_up(model, t):
    #     return model.w_dot_Delta_up[t]           >= model.w_dot[t] - model.w_dot[t-1]
    # model.constr_wtrack_up                       =  pyomo.Constraint(model.T-[1], rule=constr_wtrack_up)
    # def constr_wtrack_dn(model, t):
    #     return model.w_dot_Delta_dn[t]           >= model.w_dot[t-1] - model.w_dot[t]
    # model.constr_wtrack_dn                       =  pyomo.Constraint(model.T-[1], rule=constr_wtrack_dn)
    
    
    # # Power ramping limits
    # def constr_cramp_up(model,t):
    #     return model.w_dot_Delta_up[t]           <= model.W_dot_ramp_max*model.Delta_t
    # model.constr_cramp_up                        =  pyomo.Constraint(model.T-[1], rule=constr_cramp_up)
    # def constr_cramp_dn(model,t):
    #     return model.w_dot_Delta_dn[t]           <= model.W_dot_ramp_max*model.Delta_t
    # model.constr_cramp_dn                        =  pyomo.Constraint(model.T-[1], rule=constr_cramp_dn)
    
    
    
    # ------- Desal subsystem -------------
    # Mass balance on the desal storage
    def constr_dsmass(model,t):
        if t>1: 
            return model.m_dh[t]                 == (model.K_dhx*model.m_dot_e[t]-model.m_dot_ds[t])*model.K_t*model.Delta_t + model.m_dh[t-1]
        else:
            return model.m_dh[t]                 == (model.K_dhx*model.m_dot_e[t]-model.m_dot_ds[t])*model.K_t*model.Delta_t + model.M_dot_dh_init
    model.constr_dsmass                          =  pyomo.Constraint(model.T, rule=constr_dsmass)
    
    # Inventory limits on the desal storage 
    def constr_dsmmax(model,t):
        return model.m_dh[t]                     <= (model.M_dm_max*model.K_t)/DeltaH_dtes
    model.constr_dsmmax                          =  pyomo.Constraint(model.T, rule=constr_dsmmax)
    def constr_dsmmin(model,t):
        return model.m_dh[t]                     >= model.M_dm_min
    model.constr_dsmmin                          =  pyomo.Constraint(model.T, rule=constr_dsmmin)
    
    
    # Mass flow conversion from desal storage to distillate production 
    def constr_distillate(model,t):
        return model.v_dot[t]                    == model.K_d * model.m_dot_ds[t]
    model.constr_distillate                      =  pyomo.Constraint(model.T, rule=constr_distillate)
    
    # Range limits on desal system
    def constr_vbounds_up(model, t):
        return model.v_dot[t]                    <= model.V_dot_max
    model.constr_vbounds_up                      =  pyomo.Constraint(model.T, rule=constr_vbounds_up)
    
    def constr_vbounds_dn(model, t):
        return model.v_dot[t]                    >= model.V_dot_min
    model.constr_vbounds_dn                      =  pyomo.Constraint(model.T, rule=constr_vbounds_dn)
    
    
    # Distillate ramping tracking
    def constr_vtrack_up(model, t):
        return model.v_dot_Delta_up[t]           >= model.v_dot[t] - model.v_dot[t-1]
    model.constr_vtrack_up                       =  pyomo.Constraint(model.T-[1], rule=constr_vtrack_up)
    
    def constr_vtrack_dn(model, t):
        return model.v_dot_Delta_dn[t]           >= model.v_dot[t-1] - model.v_dot[t]
    model.constr_vtrack_dn                       =  pyomo.Constraint(model.T-[1], rule=constr_vtrack_dn)
    

# ------------------

    # ------------------
    
    
    solver = pyomo.SolverFactory('gurobi')
    results = solver.solve(model, keepfiles = True, logfile = 'solve.log')
    
    outlabs = ['t', 'C_ds', 'W_dot_gen',  'm_ch', 'm_dh', 'w_dot', 'v_dot', 'm_dot_cs','m_dot_hpt','m_dot_e','m_dot_ds', 'M_cm_max', 'M_dm_max', 'V_dot_max','C_desal']
    data_out = np.zeros((8760,15))
    for t in model.T:
        outs = [
            t, 
            model.C_ds(),
            W_dot_gen,
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
            model.C_desal()
            ]
    
        data_out[t-1] = np.array(outs)
    df_out = pd.DataFrame(data_out, columns=outlabs)
    title_str = 'dtes' + str(10*counter_dtes) + '.csv'
    df_out.to_csv(title_str, header=True)
    print(title_str)
    counter_dtes += 1


