.. _dualindnuctesscripts:

DualIndirectTES Scripts
#########################

Here is a list of relevant scripts for Model 1b and 2b.

In `neup-ies/simulations/scripts/debugging__ModelB_scripts <https://github.com/uw-esolab/neup-ies/tree/master/simulations/indtes/debugging__ModelB_scripts>`_ you can find a multitude of scripts used to debug Model 1b and 2b. These scripts are not as polished as the ones for Model 1a.


Debugging Scripts
===================

Most of these scripts are found in the ``indtes`` branch of the ``neup-ies`` repository. They will be pushed to master soon.

1. `run_PySAM_model1b.py <https://github.com/uw-esolab/neup-ies/blob/indtes/simulations/scripts/debugging__ModelB_scripts/run_PySAM_model1b.py>`_ 
	- Used to run PySAM cases for specified PC size and TES size
	- Isn't used specifically for debugging, but saves Pyomo Dispatch targets that we can then use with PySSC to actually mixed-mode debug

2. `run_PySAM_model2b.py <https://github.com/uw-esolab/neup-ies/blob/indtes/simulations/scripts/debugging__ModelB_scripts/run_PySAM_model2b.py>`_ 
	- Used to run PySAM cases for specified PC size and TES size
	- Isn't used specifically for debugging, but saves Pyomo Dispatch targets that we can then use with PySSC to actually mixed-mode debug

3. `run_PySSC_model1b.py <https://github.com/uw-esolab/neup-ies/blob/indtes/simulations/scripts/debugging__ModelB_scripts/run_PySSC_model1b.py>`_ 
	- Actually used for debugging

4. `run_PySSC_model2b.py <https://github.com/uw-esolab/neup-ies/blob/indtes/simulations/scripts/debugging__ModelB_scripts/run_PySSC_model2b.py>`_ 
	- Actually used for debugging
