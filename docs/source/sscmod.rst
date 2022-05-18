.. sscmod:

SSC Modifications
######################

For more information on the structure of a compute module in ``SSC``, see the `NREL wiki here <https://github.com/NREL/ssc/wiki/SSC-Compute-Modules#compute-module-structure>`_. Modifications to the ``SSC`` project are found `in this forked repository <https://github.com/gjsoto/ssc>`_.

Adding a new module
-------------------

By new module, I mean adding a new ``cmod`` file into the ``$DEVIR/ssc/ssc`` subfolder. There is a handy guide in the `NREL github wiki <https://github.com/NREL/ssc/wiki/SSC-Compute-Modules>`_, I have summarized some of the major steps. I have also added some linux-specific notes in addition to what is already in that wiki:

* Create a new file named ``cmod_compute_module_name.cpp`` in the ``$DEVDIR/ssc/ssc`` subfolder of the ``ssc`` project (you can do this by copying, pasting, and renaming another compute module that you wish to emulate).
* Add ``cmod_name.cpp`` to ``$DEVDIR/ssc/ssc/CMakeLists.txt``
* Compute modules must be listed twice in ``$DEVDIR/ssc/ssc/sscapi.cpp``, once in each of these two methods:

	.. code-block:: c++
	
		extern module_entry_info
		static module_entry_info *module_table[]

* Add ``cmod_name.o`` to ``$DEVDIR/ssc/build_android/Makefile-ssc`` under *Objects*
* Add ``cmod_name.o`` to ``$DEVDIR/ssc/build_ios/Makefile-ssc`` under *Objects*

.. note::
	The ``Makefile-ssc`` file was supposed to be in the ``$DEVDIR/ssc/build_linux`` subfolder but instead found it both in the ``ios`` and ``android`` subfolders. This is reflected in the steps above.
	
Adding a new solver or other file
---------------------------------

These solvers and other files would be found in the ``$DEVIR/ssc/tcs`` subfolder. They can be solvers or util files with helper methods, etc. For instance, this is where I created a new child class of the ``pt_receiver`` that acts like a simplified nuclear reactor (found at ``$DEVIR/ssc/tcs/csp_solver_nuclear.cpp``). 

The same guidance applies from the wiki: start with an existing template and start modifying it from there. Make sure to find the parent/abstract class to understand the existing member functions and constructors, then overload them correctly. 

When creating a new solver, the nomenclature seems to be to name it ``csp_solver_name`` then create both a ``.cpp`` and ``.h`` file with that same name. The header files declares member functions, variables, and structs while the ``.cpp`` file implements them. 

After creating the two corresponding files, add them to SSC by following these steps:

* Add ``csp_solver_name.o`` to ``$DEVIR/ssc/build_android/Makefile-tcs`` under *Objects*
* Add ``csp_solver_name.o`` to ``$DEVIR/ssc/build_ios/Makefile-tcs`` under *Objects*
* Add ``csp_solver_name.cpp`` and ``csp_solver_name.h`` to ``$DEVIR/ssc/tcs/CMakeLists.txt``

No additional steps are necessary (you only add ``cmod`` module files to ``ssc_api``, not solvers). 


Creating new operating mode 
----------------------------

1. Add a new class to either ``tcs/csp_solver_core.h`` or ``tcs/csp_dual_solver_core.h``
   * Should be a child of the `C_operating_mode_core` class
   * Needs a constructor (e.g., ``C_CR_OFF__PC_SB__TES_OFF__NUC_ON())``) based on base class constructor.
   * Constructor requires operating modes for CSP and PC as well as solver modes for timestep and m_dot_tes
   * Might need either `check_system_limits` or `handle_solver_error` if the base class methods need to be overloaded
  
2. Add implementations of constructors and methods in ``tcs/csp_solver_core.cpp`` or ``tcs/csp_dual_solver_core.cpp``
   
3. Add new operating mode to ``E_operating_modes`` or ``E_dual_operating_modes`` in the header file
   
4. Add new operating mode class to the ``switch-case`` in ``get_pointer_to_op_mode``
