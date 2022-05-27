.. _projectsetup:

Project Set-Up
###################################

Here is a brief overview of steps to set up the full project. 
Steps have been reproduced successfully in *Ubuntu Focal 20.04.2 LTS*.

Ultimately, through this process you will clone and build the ``SAM`` repositories: 

* `LK  <https://github.com/NREL/lk>`_  
* `WEX <https://github.com/NREL/wex>`_
* `SSC <https://github.com/gjsoto/ssc>`_
* `SAM <https://github.com/NREL/sam>`_

You will also clone:

* `googletest <https://github.com/google/googletest>`_
* `wxWidgets <https://www.wxwidgets.org/>`_

But don't worry, I have added bash scripts to this **neup-ies** repository to automate the building process.
 
.. note::
    The **neup-ies** project must be cloned into a directory alongside the other ``SAM`` projects. 
    
    It is recommended that an environment variable ``$DEVDIR`` is created to represent this parent directory.

Add repository to your PYTHONPATH
---------------------------------

First, open your *bashrc* file through an open terminal (or manually if you wish)

.. code-block:: bash

    gedit $HOME/.bashrc

then add the lines

.. code-block:: bash

    export DEVDIR=<full-path-to-parent-directory>
    export PYTHONPATH=$PYTHONPATH:$DEVDIR/neup-ies/simulations
    export git=git@github.com
    
Your ``$DEVDIR`` variable could be, for example: ``home/user_name/Documents/NE2``. Then, the ``neup-ies`` as well as the ``SAM`` directories would be placed in the shown ``NE2`` parent directory.
The github line just facilitates cloning with an SSH key.

Next, clone this repository in the parent directory. Note that you must have an `SSH key set up <https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh>`_.

.. code-block:: bash

    cd $DEVDIR
    git clone git@github.com:uw-esolab/neup-ies

    
Build Debug Version of SAM
---------------------------

* If running for the first time:
* 
	0. Make sure you delete the following directories if they exist:
		* ``$DEVDIR/build_debug``
		* ``$DEVDIR/ssc``
		* ``$DEVDIR/sam``

	1. create a CMakeLists file 
	
		.. code-block:: bash

		    cd $DEVDIR
		    touch CMakeLists.txt
		    gedit CMakeLists.txt
	
	2. Fill in the file with:
	
		.. code-block:: bash

		    cmake_minimum_required(VERSION 3.12)

		    Project(system_advisor_model VERSION 1.0.0)

		    add_subdirectory(lk)
		    add_subdirectory(wex)
		    add_subdirectory(ssc)
		    add_subdirectory(sam)

		    option(SAM_SKIP_TOOLS "Skips the sdktool and tcsconsole builds" ON)
		    option(SAM_SKIP_TESTS "Skips building tests" ON)
		    option(SAMAPI_EXPORT "Export of ssc binaries to the SAM_api directory; for Unix, compile ssc libraries for SAM_api" ON)
		    option(SAM_SKIP_AUTOGEN "Re-generate the SAMAPI files from export_config" OFF)
	
	3. Follow these particular Linux Build instructions found at `this link <https://github.com/NREL/SAM/wiki/Linux-Build-Instructions>`_:
		* Platform Requirements
		* Step 1.2 only
		* Step 2 (all steps)
			* In step 2.4, the complete path I used was ``$DEVDIR/wxWidgets-3.1.1/lib/wx-3.1.1``
	
	4. Run the bash script to build a *debug* version of ``SAM``

		.. code-block:: bash

		    cd $DEVDIR/neup-ies
		    source ./build_debug_SAM
		    
	   Bash script was created using steps found `here for mixed-mode debugging <https://github.com/uw-esolab/docs/blob/main/sam/debugSSCwithPySSC_Linux_CodeLiteIDE.md>`_.
	
	5. Verify that the project and libraries are built in the correct places:
	
		* There should be a new directory in ``$DEVDIR/build_debug``. 
		* There should be individual subdirectories for each of ``lk``, ``wex``, ``ssc``, and ``sam``
		* Check that ``$DEVDIR/build_debug/ssc/ssc/libsscd.so`` library exists
		* Check that ``googletest`` created its libraries at ``$DEVDIR/googletest/build_debug/lib``. These should be called ``libgtestd.a`` among others.
		* A CodeLite IDE workspace is created at ``$DEVDIR/build_debug/system_advisor_model.workspace``
	    
* If rebuilding a new *debug* version **OR** you already have an *export* version installed:

	0. Make sure you delete the following directories:
		* ``$DEVDIR/build_debug``
		* ``$DEVDIR/ssc``
		* ``$DEVDIR/sam``

	1. Note that the bash script at ``$DEVDIR/neup-ies/build_debug_SAM`` checks out specific branches of the ``SSC`` and ``SAM`` repositories. 
        * The script should check out a specified, stable tag of ``SAM``, if it doesn't work you could try to contact someone at NREL.
	    * The script defaults to a stable branch of my forked ``SSC`` repository, but the bash script call takes in an extra argument to override
  
	2. Run the bash script to build a *debug* version of ``SAM``

		.. code-block:: bash

		    cd $DEVDIR/neup-ies
		    source ./build_debug_SAM <optional-SSC-branch-name>

	   If you want to specify the ``SSC`` branch to check out, add an extra argument as shown with the branch name, otherwise leave that blank.

Build Export Version of SAM linked through PySAM
-------------------------------------------------

* If running for the first time:
	
	0. Make sure you delete the following directories if they exist:
		* ``$DEVDIR/build_sam_export``
		* ``$DEVDIR/build_ssc_export``
		* ``$DEVDIR/ssc``
		* ``$DEVDIR/sam``
		* ``$DEVDIR/pysam``

	1. run steps 1, 2 and 3 from the above debug section
	
	2. Run the bash script to build an *export* version of ``SAM`` and dedicated ``PySAM`` libraries

		.. code-block:: bash

		    cd $DEVDIR/neup-ies
		    source ./build_pysam

          Bash script was created using steps found `here for building PySAM with modified SSC modules <https://github.com/uw-esolab/docs/blob/main/sam/building_PySAM_using_modified_SSC.md>`_.

	3. Verify that the project and libraries are built in the correct places:
	
		* There should be a new directory: ``$DEVDIR/build_ssc_export``. 
		* There should be a new directory: ``$DEVDIR/build_sam_export``. 
		* In each of the individual subdirectories of ``lk``, ``wex``, ``ssc``, and ``sam`` there should be a ``build`` subdirectory with a ``_.a`` library
		* Check that ``$DEVDIR/build_ssc_export/ssc/libssc.so`` library exists
		* Check that ``googletest`` created its libraries at ``$DEVDIR/googletest/build/lib``. These should be called ``libgtest.a`` among others. Note this is a separate directory from the debug version
		* A CodeLite IDE workspace is created at ``$DEVDIR/build_ssc_export/sam_simulation_core.workspace``. Note that this is hardly used because currently there is no mixed-mode debugging through ``PySAM``
		* There should be .whl and .egg files in the ``$DEVDIR/pysam/dist`` directory
		* Check that ``$DEVDIR/pysam/files/libssc.so`` and ``$DEVDIR/pysam/files/libSAM_api.so`` library exists

* If rebuilding a new *export* version **OR** you already have a *debug* version installed:

	0. Make sure you delete the following directories:
		* ``$DEVDIR/build_sam_export``
		* ``$DEVDIR/build_ssc_export``
		* ``$DEVDIR/ssc``
		* ``$DEVDIR/sam``
		* ``$DEVDIR/pysam``

	1. Note that the bash script at ``$DEVDIR/neup-ies/build_pysam`` checks out specific branches of the ``pysam``, ``SSC`` and ``SAM`` repositories.  
	    * The script should check out a stable tag of ``SAM``, if it doesn't work you could try to contact someone at NREL.
	    * The script defaults to a specified, stable branch of my forked ``SSC`` repository, but the bash script call takes in an extra argument to override
	    * The script currently checks out the default ``pysam`` branch
  
	2. Run the bash script to build ``PySAM`` and an *export* version of ``SAM``

		.. code-block:: bash

		    cd $DEVDIR/neup-ies
		    source ./build_pysam <optional-SSC-branch-name>

	   If you want to specify the ``SSC`` branch to check out, add an extra argument as shown with the branch name, otherwise leave that blank.



