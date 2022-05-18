.. _dualplanttesscripts:

DualPlantTES Scripts
######################

Here is a list of relevant scripts for Model 2a (LFR + CSP + TES). 

In `neup-ies/simulations/scripts/debugging__Model2_scripts <https://github.com/uw-esolab/neup-ies/tree/master/simulations/scripts/debugging__Model2_scripts>`_ you can find a multitude of scripts used to debug Model 2a. These scripts are not as polished as the ones for Model 1a.


Debugging
===========

1. `run_PySAM_case.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/debugging__Model2_scripts/run_PySAM_case.py>`_ 
	- Used to run PySAM cases for specified PC size and TES size
	- Isn't used specifically for debugging, but saves Pyomo Dispatch targets that we can then use with PySSC to actually mixed-mode debug

2. `run_PySSCWrap_DispatchTargets.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/debugging__Model2_scripts/run_PySSCWrap_DispatchTargets.py>`_ 
	- Actually used to enter debugging mode in CodeLite or other C++ IDE


Some Data Generation
=======================
	
3. `gen__PySAM_Pyomo_sweep_for_failures__model2.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/debugging__Model2_scripts/gen__PySAM_Pyomo_sweep_for_failures__model2.py>`_ 
	- Supposed to be a parameter sweep where the static params are the user-defined table (UDT) and tariff price schedule (TPS)
	- Variable parameters we sweep over are the PC size, TES size
	- User specifies:
		- JSON script
		- tshours (as array)
		- P_ref (as array)

4. `generate_two_peak_dispatch.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/debugging__Model2_scripts/generate_two_peak_dispatch.py>`_ 
	- Creates new tariff rate schedule with two peaks (one in morning, one in late afternoon) as shown below
	.. image:: _static/dual_peak_tariffs.png
  		 :target: _static/dual_peak_tariffs.png

5. `generate__load_profiles.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/debugging__Model2_scripts/generate__load_profiles.py>`_ 
	- Generates a full year of simulation profiles for model2a
	- User specifies:
		- JSON script
		- tshours (single value)
		- P_ref (single value)

	        
Post Processing / Plotting
============================

6. `plot__load_profiles.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/debugging__Model2_scripts/plot__load_profiles.py>`_ 
	- Creates the violin plots for Twin Peak Results
	- Splits distributions between Winter vs Summer, Weekday vs Weekend
	- Calls on a util class called ``NuclearTESLoadProfiles`` which hosts a lot of methods used to make the violin plots (some variations to deal with CSP profile)
	- Example plot:
	.. image:: _static/model2a_loadprof.png
	        :target: _static/model2a_loadprof.png	

7. `plot__PySAM_Pyomo_sweep_for_failures__model2.py <https://github.com/uw-esolab/neup-ies/blob/master/simulations/scripts/debugging__Model2_scripts/plot__PySAM_Pyomo_sweep_for_failures__model2.py>`_ 
	- plotting results from parameter sweep files
	- Param sweep files are found in the ``neup-ies/simulations/outputs`` subdirectory
		- usually called something like "failureModes_PySAM__model2_Hamilton_560_tariffx1__2022_02__pyomo_1__horizon_24_48__TES_[2,14]__PC_[600,850]__CSP_400"



