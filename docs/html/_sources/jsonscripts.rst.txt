.. _jsonscripts:

JSON Scripts
#################################################

JSON scripts in `neup-ies/simulations/json <https://github.com/uw-esolab/neup-ies/tree/master/simulations/json>`_ are defined as nested "dictionaries" in the following format:

.. image:: _static/json_structure.png
   :target: _static/json_structure.png
   
PySAM_inputs
=============

These inputs are typically either parameters useful in the Pyomo optimization that are not used within SSC or filepaths to data files that we want to extract data and set as vector/array inputs to SSC.


PySAM_outputs
==============

These are keywords that are used within the NE2 modules to generate an .csv output file (optional).


SSC_inputs
============

As of Python 3.6, dictionary keys preserve order when retrieved. Therefore the ordering of entries in the `SSC_inputs` dictionary needs to be in a specific format.

We assume that there are data entries for each of three modules:
- a `Plant` module (`tcsmolten_salt`, `nuclear_tes`, etc.)
- a `Grid` module (typically `grid`)
- a `Financial` module (typically `single_owner`)

The **first** entry of `SSC_inputs` should be:

.. code-block::

	"compute_module_0" : <SSC_module_name>, 

for whatever module you are trying to simulate. All entries after that belong to the `Plant` module. 

The **next** entry (after the last `Plant` data entry is specified) should be:

.. code-block::

	"compute_module_1" : <SSC_grid_name>, 

followed by `Grid`-specific data entries.

Repeat the syntax for `Financial` module data entries.

**Finally**, the last entry of the `SSC_inputs` dictionary in the JSON script should be:

.. code-block::

	"number_compute_modules" : <N> 
	}
	
where *N* is the number of compute modules (typically 3).
This syntax ensures that we can run the `PySSCWrapper` correctly for debugging.


