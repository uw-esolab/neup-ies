.. _dualindnuctesstatus:

Current Status of Models 1b and 2b
###################################

Outlined here is the status of the Model 1b and 2b code framework.


Status of Model Components
================================

- **Dispatch Model**: [*complete*] a Pyomo dispatch MILP has been created and runs correctly within the project. It sets up an MILP and solves it when called.

- **SSC Model**: [*in progress*] new classes have been added to SSC suite and work but calling through PySAM is currently not working (only through PySSC). ``mdot_solver`` also needs to be overloaded with correct stream calculations.

- **NE2 Model**: [*complete*] a Python interface class has been created and runs correctly within the project. It creates dispatch model instances, solves them, converts results to dispatch targets, and calls on SSC to execute at every timestep of the time horizon. 


Framework within SSC
================================

1. **Compute module (cmod)**: *cmod_nuclear_mspt_indirect_tes*

	- This part works
	- major modifications include new parameters and setting them correctly within member objects
	
2. **csp_solver**: *csp_dual_indirect_solver_core*

	- This part mostly works
	- needs some debugging within the estimation of plant component modes and operating mode selection

3. **csp_solver_mono_eq**: *csp_dual_indirect_solver_mono_eq_methods*

	- This part is not finished
	- Need to overload mass flow solver - calculate the correct streams and temperatures and set the proper inputs within the PC and TES object (mix the streams, do HX calculations, etc.)

3. **Plant Component solvers**: within the compute module, an object instance of all plant components are initialized with corresponding SSC inputs. Then these objects are used as inputs for the ``csp_solver``.

	- TES: *csp_solver_two_tank_tes*
		- No changes
	- PC: *csp_solver_pc_Rankine_indirect_224*
		- Some changes enacted, debugging
	- CSP: *csp_solver_mspt_collector_receiver*
		- Receiver: *csp_solver_mspt_receiver_222*
		- Heliostatfield: *csp_solver_pt_heliostatfield*
			- No changes
	- LFR: *csp_solver_nuclear_plant*
		- Nuclear: *csp_solver_nuclear*
			- Some changes enacted, debugging
		
