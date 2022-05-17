.. _sscoverview:

SSC Overview
##############

There are a lot of levels to SSC but the major ones we highlight are:

1. **Compute module (cmod)**: these are the staging areas for gathering inputs, initializing different C++ classes (TES, PC, CSP, etc.)

	- *cmod_tcsmolten_salt*: used for CSP + TES
	- *cmod_nuclear_tes*: used for **Model 1a** (LFR + TES)
	- *cmod_nuclear_mspt_tes*: used for **Model 2a** (LFR + CSP + TES direct)
	- *cmod_nuclear_mspt_indirect_tes*: used for **Model 1b, 2b** (LFR + CSP + TES indirect)
	
2. **csp_solver**: this sets up methods and member attribute objects used to carry out engineering model calculations

	- *csp_solver_core*
	- *csp_dual_solver_core*
	- *csp_dual_indirect_solver_core*

3. **csp_solver_mono_eq**: here are extra methods used to solve mass flow, defocus, etc. at each timestep. Also calls individual plant component solvers.

	- *csp_solver_mono_eq_methods*
	- *csp_dual_solver_mono_eq_methods*
	- *csp_dual_indirect_solver_mono_eq_methods*

3. **Plant Component solvers**: within the compute module, an object instance of all plant components are initialized with corresponding SSC inputs. Then these objects are used as inputs for the ``csp_solver``.

	- TES: *csp_solver_two_tank_tes*
	- PC: *csp_solver_pc_Rankine_indirect_224*
	- CSP: *csp_solver_mspt_collector_receiver*
		- Receiver: *csp_solver_mspt_receiver_222*
		- Heliostatfield: *csp_solver_pt_heliostatfield*
	- LFR: *csp_solver_nuclear_plant*
		- Nuclear: *csp_solver_nuclear*
		
	.. note::
		
		The LFR is set to mimic the CSP. Within the ``csp_solver`` a member object is created for the ``collector_receiver`` which itself has two members: a receiver and a heliostatfield object.
		
		The LFR also has a two-tiered system: a ``nuclear_plant`` object which has a member ``nuclear`` object. 

For a more detailed view as to which csp solvers we use with each Model, see :ref:`here <modulecombinations>`.
