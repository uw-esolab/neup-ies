.. _runningpysam:

Running PySAM
#####################

There are two ways to run simulations within this project: using PySAM (**recommended**) and using :ref:`PySSC <runningpyssc>` (**ideally for debugging**). 

The main way to run a full simulation in SSC **with dispatch optimization** is to use the NE2 modules.
These are found in this repository under `neup-ies/simulations/modules <https://github.com/uw-esolab/neup-ies/tree/master/simulations/modules>`_. 

A number of example scripts for running a simulation can be found under `neup-ies/simulations/scripts <https://github.com/uw-esolab/neup-ies/tree/master/simulations/scripts>`_. There should be a guide on module-specific scripts within this documentation.

Here is a general outline of how to run a simulation.

0. **Build SSC** in *export* mode. Instructions found :ref:`here <projectsetup>`:. 

1. **PYTHONPATH**: Make sure the ``neup-ies/simulations/`` directory is in your Python path, either in your terminal or whatever IDE you are using. This is slightly automated in current scripts by calling:

	.. code-block:: bash

	    import sys
	    sys.path.append('..')

2. **Select an NE2 module** you want to run. For example, if you want to run an LFR + TES model (**Model 1a**) you would select the ``NuclearTES`` module. Import the module at the top of your script:

	.. code-block:: python

	    import modules.NuclearTES as NuclearTES

3. **Initialize** an object instance of the module class as follows:

	.. code-block:: python

	    nuctes = NuclearTES.NuclearTES()

   Note that there are several input parameters you can specify, though each module has default values preset. 
   
        - plant_name (str) – name of SSC module to run

        - json_name (str) – name of JSON script with input data for module (just the name, don't specify directory)
        
        	.. note::
	   
		   JSON scripts can be found `neup-ies/simulations/json <https://github.com/uw-esolab/neup-ies/tree/master/simulations/json>`_. 
		   For more information on the structure of the JSON scripts, see :ref:`this guide <jsonscripts>`:.
	   

        - is_dispatch (bool) – boolean, if True runs Pyomo dispatch optimization

        - dispatch_time_step (int) – time step for dispatch (hours)

        - log_dispatch_targets (bool) – boolean, if True logs dispatch targets calculated by Pyomo at each segment
        	
        	.. note::
	   
	   	   The ``log_dispatch_targets`` input, if set to True, saves all dispatch target time series that are sent to SSC from the Pyomo optimization model.
	   	   These time series are saved to unique files ending in `.dispatchtargets` whose filenames are hash strings created from all JSON script parameter inputs.
	   	   These files are found in `neup-ies/simulations/outputs <https://github.com/uw-esolab/neup-ies/tree/master/simulations/outputs>`_. 

        - exec_debug (bool) – boolean, allows execution in “debug” mode that times out exec call

        - exec_timeout (float) – amount of time in seconds to timeout an exec call

	       .. note::
	   
	   	   ``exec_debug`` and ``exec_timeout`` were intended to help debug sessions were calls to the root SSC program would get stuck.
	   	   These inputs have limited success in helping debug those sessions...
	   	   
4. **Optional Changes for Sweeps**: After Step 3, the JSON script parameters should be loaded into the object. If you want to make any alterations (say you are running a loop or parameter sweep), there are some extra lines you should run before runnning the simulation to update all stored dictionaries.

	- If the change is to a ``PySAM_dict`` parameter (say, the ``ssc_horizon`` or ``pyomo_horizon`` time):

		.. code-block:: python

		    # horizons
		    nuctes.PySAM_dict['ssc_horizon']   = sscH
		    nuctes.ssc_horizon   = sscH * nuctes.u.hr
		    nuctes.PySAM_dict['pyomo_horizon'] = pyoH
		    nuctes.pyomo_horizon = pyoH * nuctes.u.hr

		    # saving/updating PYSAM dict to nuctes
		    nuctes.dispatch_wrap = nuctes.create_dispatch_wrapper( nuctes.PySAM_dict )
		    
	  The last line will re-create a new instance of a the Dispatch Wrapper class that communicates with the Pyomo model.
	
	- If the change is to an ``SSC_dict`` parameter (say, the design point `P_ref` or `tshours`):
	
		.. code-block:: python

		    nuctes.SSC_dict['P_ref'] = Pref
		    nuctes.SSC_dict['tshours'] = tshours
		    nuctes.dispatch_wrap.set_design()
		    
	  The last line updates the existing Dispatch Wrapper class (which is a member of the module class) with the new design point. 
		    
	
5. **Run Simulation**: 

	.. code-block:: python

	    nuctes.run_sim(  )
	    
	Here are some of the input parameters to this method:


        - run_loop (bool) – if true, runs simulation in segments. else, runs simulation all at once
        
        	.. note::
	   
	   	   If running with dispatch optimization set to True, make sure ``run_loop`` is set to True.

        - export (bool) – if true, exports results to an Excel sheet

        - filename (str) – name for Excel sheet saved to the /outputs directory

        - overwrite_dispatch_targets (bool) – if true, overwrites the current stored dispatch target file
        
                .. note::
	   
	   	   Normally, if ``log_dispatch_targets`` is set to True in the init of the object, the object checks to see if a dispatch targets file 
	   	   exists for the unique hash signature created out of the JSON script input parameters.
	   	   
	   	   If the file exists, it doesn't save the current simulations's results to that file.
	   	   
	   	   If ``overwrite_dispatch_targets`` is set to True here, it overrides the previous functionality and overwrites the existing file.




