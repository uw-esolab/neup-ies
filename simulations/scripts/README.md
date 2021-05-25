# Scripts for Simulations

This folder is for Python simulaton scripts.

Because this can get clogged up fast, for now I just added all *.py files in this folder to .gitignore with 3 exceptions. 

If you'd like to upload a specific Python file, override the .gitignore by:
>git push -f

TODO: come up with a naming convention?

## Running PySAM scripts
Be sure to follow the instructions in the directory above this one to build PySAM correctly. 

Use `PySAM_model1.py` as a template to run PySAM.

## Debugging SSC
Currently we have no way of debugging through PySAM.
Instead, we can run things through PySSC and debug the C++ code from there. 

- In the **PySSC_scripts** subdirectory, the `PySSC.py` script contains a class to class `SSC` C++ code directly and some helper methods.
- The **../util/PySSCWrapper.py** script contains a class to wrap PySSC and reads/assigns parameters from a user-selected JSON script similar to how PySAM works. 
- Use `PySSC_scripts/run_PySSCWrap_model1.py` as a template for debugging SSC. Read [the following docs, especially Section 5](https://github.com/uw-esolab/docs/blob/main/sam/debugSSCwithPySSC_Linux_CodeLiteIDE.md) for more detailed instructions on mixed-mode debugging (Python calling C++). 

