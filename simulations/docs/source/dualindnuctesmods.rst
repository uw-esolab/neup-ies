.. _dualindnuctesmods:

Modifications within Model 1b and 2b
######################################

Outlined here are some major modifications of the Model 2a SSC code.

csp_solver Changes
=====================

Major changes to ``csp_solver`` include:

	- creating a new ``csp_dual_indirect_solver`` that is a child of ``csp_dual_solver``

	- Added new input dispatch target: fraction of nuclear output going to TES
	
	- New input parameters:
	
		- Steam-to-salt HX effectiveness
		
		- Salt-to-steam HX effectiveness
	
	- Two HTF fluids:
		
		- CSP and TES use salt as before
		
		- PC and LFR use #4 - steam
	
	- Calculating design points using enthalpies rather than specific heat formula
	

csp_solver_pc_Rankine_indirect_224 Changes
==============================================

	- Calculating design points using enthalpies rather than specific heat formula
	
	- Using ``water_props`` class instead of ``HTFProperties`` to get steam properties


csp_solver_nuclear Changes
==============================================

	- Calculating design points using enthalpies rather than specific heat formula
	
	- Using ``water_props`` class instead of ``HTFProperties`` to get steam properties


csp_solver_mono_eq_methods Changes
====================================

Need to:

	- Add 2 HX models within mass flow solvers to calculate:
		
		- salt mass flow from LFR to hot tank
		
		- steam mass flow through TES dispatch
	
	- Recombine steam mass flows through LFR and through TES before calling PC method 
	
	
