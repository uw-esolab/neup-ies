# Scripts for Simulations

This folder is for Python simulaton scripts.

Because this can get clogged up fast, for now I just added all *.py files in this folder to .gitignore with 3 exceptions. 

If you'd like to upload a specific Python file, override the .gitignore by:
>git push -f

TODO: come up with a naming convention?

## Running PySAM scripts
Be sure to follow the instructions in the **neup-ies** main directory to build PySAM correctly. 
- All PySAM modules are defined in the **neup-ies/simulations/modules** subdirectory
- Use `neup-ies/simulations/scripts/PySAM_model1.py` as a template to run PySAM.
- The corresponding CodeLite workspace for the C++ SSC project is found in `../build_ssc_export/sam_simulation_core.workspace` (above the **neup-ies** directory) 

## Debugging SSC
Currently we have no way of debugging through PySAM.
Instead, we can run things through PySSC and debug the C++ code from there. 
Be sure to follow the instructions in the **neup-ies** main directory to make sure a debugging version of all projects is built alongside the export ones for PySAM. 

- In the **PySSC_scripts** subdirectory, the `PySSC.py` script contains a class to class `SSC` C++ code directly and some helper methods.
- The `neup-ies/util/PySSCWrapper.py` script contains a class to wrap PySSC and reads/assigns parameters from a user-selected JSON script similar to how PySAM works. 
- Use `neup-ies/simulations/scripts/PySSC_scripts/run_PySSCWrap_model1.py` as a template for debugging SSC. 
- Read [the following docs, especially Section 5](https://github.com/uw-esolab/docs/blob/main/sam/debugSSCwithPySSC_Linux_CodeLiteIDE.md) for more detailed instructions on mixed-mode debugging (Python calling C++). 
- The corresponding CodeLite workspace for the C++ SSC project is found in `../build_debug/system_advisor_model.workspace` (above the **neup-ies** directory) to use when debugging

