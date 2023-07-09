#######################################
Contents:

The Aryan-python-folder contains:

1.) rankine-on-design-ttd-06-24-23.py
	- The most recent model for simulating supercritical steam rankine cycle in on-design state. 
	- Does not include any solar loop integration.
	- Based off of \Aryan-python-model\white-supplement.zip\models\supercritical-steam-rankine\rankine-on-design.EES

2.) components_v2.py
	- component model library used in conjunction with on-design model

3.) state.py
	- python library used to fix state

4.) golden section search folder
	- folder of files containing model based on old golden section search.

5.) white-supplement
	- EES models of Brian White, used as base for all python models

#######################################
Rankine cycle model:

rankine-on-design-ttd-06-24-23.py
	- OFF design implemented 
	- TTD and DT values are set to acheive scaled UA values
	- Framework setup to include 100% charging and discharging --> Need to add solar loop states and change equation pipeline.
	- Implicit model. Convergence based currently on m[15] over cycles. Some issues seen.

Future steps/goals
	- Add convergence status based on CoE --> model iterates till CoE is balanced
	- Turn entire code into a modular solver function to allow for parametric studies
	- Switch solver to include a cost function and separate optimizer routing that uses scipy's inbuilt optimizer instead of manually programming optimization algorithm.
