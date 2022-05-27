.. _nuctesscripts:

NuclearTES Scripts
######################

Here is a list of relevant scripts for Model 1a (LFR + TES). A lot of these scripts were used to generate data for the Energies paper `(link here) <https://doi.org/10.3390/en15103599>`_.

In `neup-ies/simulations/scripts/output_scripts <https://github.com/uw-esolab/neup-ies/tree/master/simulations/scripts/output_scripts>`_ you can find a multitude of scripts to both generate data and convert existing data into figures.

Within that folder, there are some subfolders:

 - **ANS2021Winter**: scripts used to generate data for the ANS 2021 Winter Meeting back in November 2021
 
 - **Model1Sweeps**: data generation scripts, data automatically saved to the ``simulation/outputs`` directory
 
 - **Model1Figures**: scripts to take data from ``simulation/outputs`` and generate figures


Generating Data Used in Paper
===============================

Within **Model1Sweeps**...

1. `generate_exaggerrated_tariff_schedules.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/generate_exaggerrated_tariff_schedules.py>`_ 
	- Creates a tariff schedule from SAM Tariff rate with amplified peaks and troughs. 
	- Amplification factor can be specified by user

2. `paramSweep_stat_UDT_TPS__var_PC_TES_TurbCost.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/paramSweep__stat_UDT_TPS__var_PC_TES_TurbCost.py>`_
	- Supposed to be a parameter sweep where the static params are the user-defined table (UDT) and tariff price schedule (TPS)
	- Variable parameters we sweep over are the PC size, TES size
	- This uses Cory's empirical formula for overdesign turbine cost
	- User specifies:
		- JSON script
		- tshours (as array)
		- P_ref (as array)
	- Related scripts for stacking results from multiple runs together:
		a. `stitch_results__P_ref_multiruns.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/stitch_results__P_ref_multiruns.py>`_ 
			- stacks results from multiple runs with different P_ref ranges
		b. `stitch_results__tshours_multiruns.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/stitch_results__tshours_multiruns.py>`_  
			- stacks results from multiple runs with different tshours ranges
		- Used these scripts because for a single parameter sweep, I would split the runs over several consoles/cores as a very manual parallelization

3. `paramSweep_sensitivity_studies_TurbCost.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/paramSweep_sensitivity_studies_TurbCost.py>`_
	- Supposed to be a parameter sweep where the static params are the user-defined table (UDT) and tariff price schedule (TPS)
	- Variable parameters we sweep over TES costs
	- User specifies:
		- JSON script
		- tshours (discrete set)
		- P_ref (discrete set)

Within **Model1Figures**...

4. `generate__load_profiles.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/generate__load_profiles.py>`_ 
	- Generates a full year of simulation profiles
	- User specifies:
		- JSON script
		- tshours (single value)
		- P_ref (single value)

		
Recreating Figures from Paper
==============================

All outputs from the previous section are saved to `simulations/outputs/model1a_energies_paper <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/outputs/model1a_energies_paper>`_ 
Within **Model1Figures**...

1. `plot__load_profiles.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/plot__load_profiles.py>`_ 
	- Creates the violin plots for SAM Tariff rate results
	- Splits distributions between Winter vs Summer, Weekday vs Weekend
	- Calls on a util class called ``NuclearTESLoadProfiles`` which hosts a lot of methods used to make the violin plots
	
2. `plot__load_profiles__caiso.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/plot__load_profiles__caiso.py>`_ 
	- Creates the violin plots for CAISO Tariff rate results
	- Splits distributions between Winter vs Summer, Weekday vs Weekend
	- Calls on a util class called ``NuclearTESLoadProfiles`` which hosts a lot of methods used to make the violin plots

3. `paramSweep__tariff_exaggerated.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/paramSweep__tariff_exaggerated.py>`_ 
	- Creates the relative PPA heatmap plots for SAM Tariff rates of PC size vs TES size

4. `paramSweep__caiso.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/paramSweep__caiso.py>`_ 
	- Creates the relative PPA heatmap plot for CAISO Tariff rates of PC size vs TES size

5. `paramSweep__plot_sensitivity_on_TESCost.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/paramSweep__plot_sensitivity_on_TESCost.py>`_ 
	- Creates the 3x3 plots of TES Cost sensitivity
	
6. `pricing_multipliers_CAISO.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/output_scripts/Model1Sweeps/pricing_multipliers_CAISO.py>`_ 
	- Creates the pricing multipliers as time series plot
	
	
