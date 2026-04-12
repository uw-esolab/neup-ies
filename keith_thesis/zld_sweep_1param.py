import pyomo.environ as pyo
import pandas as pd
import numpy as np
import os
import csv
import glob
from zld_dynamic import build_model



def run_1d_sensitivity():
    
    
    N_design    = 24
    N_dispatch  = 8760

    # Load Data
    eff_design    = pd.read_csv("avg_day_efficiency.csv")["Average_Efficiency"].to_numpy()[:N_design]
    elec_design   = pd.read_csv("avg_day_lmp.csv")["LMP"].to_numpy()[:N_design]
    eff_dispatch  = pd.read_csv("year_efficiency.csv")["Efficiency"].to_numpy()[:N_dispatch]
    elec_dispatch = pd.read_csv("year_lmp.csv")["LMP"].to_numpy()[:N_dispatch]

    solver = pyo.SolverFactory("gurobi")
    solver.options['NonConvex'] = 2


    # define sweep ranges
    # sweep lithium price multiplier
    sweep_name     = 'lithium'
    sweep_values   = [0.25, 0.5, 1.0, 2.5, 5.0, 7.5, 10.0, 20.0]
    
    
    # sweep water price multiplier
    # sweep_name   = 'water'
    # sweep_values = [0.25, 0.5, 1.0, 1.5, 2.0, 5.0]
    
    
    # sweep cold tes cost multiplier
    # sweep_name   = 'tes_cold'
    # sweep_values = [0.25, 0.5, 1.0, 2.5, 5.0, 10.0]
    
    
    # stretched lmp files
    price_schedule_files = sorted(
        glob.glob(os.path.join("stretched_lmp_cases", "year_lmp_top*_bottom*.csv"))
)

    
    results_dir = os.path.join("sweeps", f"sweep_1var_{sweep_name}")
    os.makedirs(results_dir, exist_ok=True)
    
    # 
    
    summary_data = []

    # loop over multipliers
    for mult in sweep_values:
        print(f"\n>>> Solving case: {sweep_name} multiplier = {mult}")
        current_prices = {sweep_name: mult}

        # design phase
        model_d = build_model(N_design, elec_design, eff_design, price_mults=current_prices)
        solver.solve(model_d, tee=False)


        # extract sizing
        caps      = {p: pyo.value(model_d.v_dot_cap[p]) for p in model_d.processes}
        hot_max   = pyo.value(model_d.m_hot_tes_max)
        cold_max  = pyo.value(model_d.m_cold_tes_max)
        y_p       = {p: pyo.value(model_d.y_process_active[p]) for p in model_d.processes}
        y_l       = {l: pyo.value(model_d.y_link_active[l]) for l in model_d.links}


        # dispatch phase
        model_disp = build_model(N_dispatch, elec_dispatch, eff_dispatch, price_mults=current_prices)
        
        
        # fix variables
        for p in model_disp.processes:
            model_disp.v_dot_cap[p].fix(caps[p])
            model_disp.y_process_active[p].fix(round(y_p[p]))
            
        for l in model_disp.links:
            model_disp.y_link_active[l].fix(round(y_l[l]))
            
        model_disp.m_hot_tes_max.fix(hot_max)
        model_disp.m_cold_tes_max.fix(cold_max)

        res_disp = solver.solve(model_disp, tee=False)


        if res_disp.solver.termination_condition == pyo.TerminationCondition.optimal:

            mult_str = str(mult)
            hourly_filename = os.path.join(
                results_dir,
                f"dispatch_{sweep_name}_{mult_str}.csv")
            
            
            with open(hourly_filename, "w", newline="") as f:
                w = csv.writer(f)

                header = ["hour", "hot_tes_kWh", "cold_tes_kWh", "elec_sold_kWe", "elec_gen_kWe", "m_dot_extract_kg/s", "li_recovered_kg"]
                for p in model_disp.processes:
                    header.append(f"v_in_{p}")
                    header.append(f"v_dil_{p}")
                w.writerow(header)

                for t in model_disp.T:
                    
                    hot_tes_kWh = (pyo.value(model_disp.m_hot_tes[t]) * model_disp.Cp_salt.value * model_disp.Delta_T_hx1_salt.value/3600)
                    cold_tes_kWh = (pyo.value(model_disp.m_cold_tes[t]) * model_disp.Delta_H_hx2_water.value/3600)
                    
                    row = [
                        t,
                        hot_tes_kWh,
                        cold_tes_kWh,
                        pyo.value(model_disp.elec_sold[t]),
                        pyo.value(model_disp.w_dot_gen[t]),
                        pyo.value(model_disp.m_dot_extract[t]),
                        pyo.value(model_disp.li_recovered_ed[t]), ]
                    
                    for p in model_disp.processes:
                        row.append(pyo.value(model_disp.v_dot_in[t, p]))
                        row.append(pyo.value(model_disp.v_dot_dil[t, p]))
                    
                    w.writerow(row)

            print(f"Saved: {hourly_filename}")


        else:
            print(f"Dispatch model failed for {sweep_name} = {mult}")


    df_summary = pd.DataFrame(summary_data)
    df_summary.to_csv(f"{results_dir}/sweep_summary_1d.csv", index=False)
    print(f"\nSweep complete. Summary saved to {results_dir}/sweep_summary_1d.csv")

if __name__ == "__main__":
    run_1d_sensitivity()




























