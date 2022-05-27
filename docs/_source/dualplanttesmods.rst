.. _dualplanttesmods:

Modifications within Model 2a
###############################

Outlined here are some major modifications of the Model 2a SSC code.

csp_solver Changes
=====================

Major changes to ``csp_solver`` include:

	- creating a new ``csp_dual_solver`` that is a child of ``csp_solver``
	
	- In the original, the constructor takes in a single ``collector_receiver`` object (in Model 1a, we use a child class called ``nuclear_plant``)
	
		- Here, we have a new constructor that takes it a ``collector_receiver`` **AND** a ``nuclear_plant`` object
		
	- A "nuclear" version of a lot of CSP input parameters are created
	
	- Added calls specifically to ``nuclear_plant`` estimation methods before operating mode selection
	
		- Some if-else statements have been added to check if targets from Pyomo are reasonable
		
		- Now accounting for mixing of LFR and CSP outlet streams to hot tank
		
	- Overloading the ``find_operating_modes`` method 
	
		- have to add option to select new operating modes
		
		- account for dual heat and mass flow from both CSP and LFR
		
		- here is a gigantic diagram of the current state of operating mode selection: 
		
		.. image:: _static/model2a_operatingmodeselect.png
    			:target: _static/model2a_operatingmodeselect.png
    			
	- Logging outputs modified to account for added LFR stream


New Operating Modes
=====================

Below is a diagram of all new operating modes (dotted border) categorized by the corresponding CSP, PC, and TES modes that already exist.

.. image:: _static/model2a_newopmodes.png
   :target: _static/model2a_newopmodes.png
   

csp_solver_mono_eq_methods Changes
====================================

Major changes to ``csp_solver_mono_eq_methods`` include:

	- creating a new ``csp_dual_solver_mono_eq_methods`` that is a child of ``csp_solver_mono_eq_methods``
	
	- Mass flow solver now accounts for flow mixing between LFR and CSP output before entering hot tank

	.. image:: _static/model2a_mdotsolver.png
   		:target: _static/model2a_mdotsolver.png
