.. _highlevelcode:

Code Overview
#################################################

State of affairs as of 17 May 2022. 

Current implementation does not include RAVEN functionality. RAVEN can be wrapped around this whole code.

.. image:: _static/NE2_full_codeflow_asOfMay17.png
   :target: _static/NE2_full_codeflow_asOfMay17.png


The different aspects of the codeflow are described as follows:

0. **JSON Scripts**: these contain input information necessary to run SSC, including path links to data files used within either SSC or the NE2 modules. They also contain extraneous parameters used within the Pyomo optimization. More information can be found in :ref:`General Project Structure - JSON Scripts <jsonscripts>`.

1. **Python Scripts**: essentially a staging area for the simulations. In some sort of script, we initialize the necessary NE2 modules and link them to a JSON script. This can also happen within a terminal or within some RAVEN interface. General info about what goes in this staging process can be found :ref:`Guides to Building/Running - Running PySAM <runningpysam>`.

2. **NE2 Modules**: within these classes, we create PySAM modules that communicate with the desired SSC module. We also loop over the SSC time horizon (1 day) for a full year: during each day we stagger calls to the pyomo dispatch optimization and then to SSC. Here we log all simulation results. 

3. **Dispatch Parameter Wrap**: here are methods that create and/or update a parameter dictionary with inputs to the dispatch optimization problem.

4. **Dispatch Class**: here are methods that actually create the Pyomo MILP concrete model (adds parameters, variables, constraints, and an objective function). Contains a method to solve the MILP model.

5. **Dispatch Outputs Class**: this class converts outputs from dispatch optimization into targets for the SSC inputs. 

6. **SSC**: this is the SAM Simulation Core engineering model. Information on the engineering model is specific to the intended IES model. 

7. **PostProcessing**: there are some classes in `neup-ies/simulations/util <https://github.com/uw-esolab/neup-ies/tree/master/simulations/util>`_ that are useful for plotting and displaying simulation outputs. This can also be replaced with some RAVEN functionality if desired. 


For a more detailed diagram of how to combine classes and modules from items **2**-**6**, see :ref:`General Project Structure - NE2, Dispatch, and SSC Module Combinations <modulecombinations>`.

