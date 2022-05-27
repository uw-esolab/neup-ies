.. ne2mod:

NE2 Modules 
#############

There exist three distinct sets of modules:

    * SSC modules   : these are the `cmod_<name>.cpp` files from the SSC repository written in C++ (e.g. cmod_tcsmolten_salt, cmod_nuclear_tes)
    * PySAM modules : these are PySAM wrappers for the SSC modules (also Grid, Financial modules, etc.) which are created from the `build_pysam` bash script
    * NE2 modules   : these are Python classes that create PySAM modules (e.g., NuclearTES, SolarTES)

Here is a general diagram of basic NE2 module structure (to be updated later):

.. image:: _static/NE2_modules_overview.png
   :target: _static/NE2_modules_overview.png
