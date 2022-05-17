.. _runningpyssc:

Running PySSC
################

There are two ways to run simulations within this project: using :ref:`PySAM <runningpysam>`. (**recommended**) and using PySSC (**ideally for debugging**). 

Instead of calling the NE2 modules, we call instead use a single class called ``PySSCWrapper`` found in `neup-ies/simulations/util <https://github.com/uw-esolab/neup-ies/tree/master/simulations/util>`_. 

.. note::

    This is **not** set up for using with Pyomo optimization. 
    You can, however, run a simulation with dispatch optimization using the NE2-PySAM script which will save a dispatch target time series.
    PySSC runs can then automatically load those dispatch target files (if the JSON inputs are the same!) and use them as a vector input. 
    
    
A number of example scripts for running a simulation can be found under `neup-ies/simulations/scripts/PySSC_scripts <https://github.com/uw-esolab/neup-ies/tree/master/simulations/scripts/PySSC_scripts>`_. Notably, `run_PySSCWrap_model1 <https://github.com/uw-esolab/neup-ies/tree/master/simulations/scripts/run_PySSCWrap_model1.py>`_ is the cleanest. 

Here is a general outline of how to run a simulation.

0. **Build SSC** in either *export* or *debug* mode. Instructions found :ref:`here <projectsetup>`:. 

1. **PYTHONPATH**: Make sure the ``neup-ies/simulations/`` directory is in your Python path, either in your terminal or whatever IDE you are using.

2. **Optional: Print PID**: To help in setting up mixed-mode debugging, you can print out the process ID number so you can attach a C++ IDE to the current process:

	.. code-block:: python
	
	   pid = os.getpid()
	   print("PID = ", pid)
	   
3. **Initialize** the PySSCWrapper object

	.. code-block:: python
	
	   pw = PySSCWrapper()
	   
   Note that there are two input parameters you can specify: 
   
        - json_name (str) – name of JSON script with input data for module (just the name, don't specify directory)
        
        	.. note::
	   
		   JSON scripts can be found `neup-ies/simulations/json <https://github.com/uw-esolab/neup-ies/tree/master/simulations/json>`_. 
		   For more information on the structure of the JSON scripts, see :ref:`this guide <jsonscripts>`:.
	   

        - is_debug (bool) – boolean, if True calls on SSC through debug .so library file, else chooses the export library .so file

4. **Run Simulation**: 

	.. code-block:: python

	    pw.run_sim()

   Note that there is one input parameter you can specify: 
   
        - run_dispatch_targets (bool) – if True, loads dispatch target time series from saved hash file whose name matches JSON inputs (only if it exists).
        

Mixed-Mode Debugging
======================
        
For a detailed guide on setting up mixed-mode debugging (between Python and C++) check `here <https://github.com/NREL/ssc/wiki/Debug-SSC-with-PySSC-Linux-CodeLiteIDE>`_. I also wrote a separate one for Eclipse IDE that is also in that wiki. 

Mixed-mode debugging would be possible in Microsoft Visual Studio Code but have not personally figured that out yet. 

