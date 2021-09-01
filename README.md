# neup-ies

## Integrated Solar & Nuclear Cogeneration of Electricity & Water using the sCO2 Cycle

This is a private repository for sharing code and project files.

*Note from Gabriel: If cloning onto a local machine, I recommend that this repository be cloned into a directory alongside the SAM projects (LK, WEX, SSC, SAM). That is, the **neup-ies** folder should be in the same directory as **lk**, **ssc**, **googletest**, etc. I refer to this directory as `$DEVDIR`.*
<br/><br/>

## Building Documentation

Documentation written using Sphinx <http://sphinx-doc.org/>_. 
Follow these steps to build the documentation for this project!

Make sure you have the proper packages installed in your Python environment:
    
    pip install Sphinx
    pip install numpydoc
    pip install sphinxcontrib-napoleon

First make sure you are in the correct conda environment. Then, make sure you run:

    cd neup-ies/simulations/docs
    sphinx-apidoc -M -f -o source/ ../../simulations/
    make html

The final html file will be located at ``neup-ies/simulations/docs/build/html/index.html``.

# Running SSC simulations in Linux
Everything in the bash scripts has been tested on Ubuntu Focal 20.04.2 LTS.

## Add repository to your PYTHONPATH
Make sure to add the simulations subdirectory to your Python path. First open your `bashrc` file, could do it from the command line:
```
gedit $HOME/.bashrc
```
and add the line
```
export PYTHONPATH=$PYTHONPATH:$DEVDIR/neup-ies/simulations
```
where `$DEVDIR` is the parent directory where `neup-ies` is located. 
<br/><br/>

## Building SAM
Before building either the debug or export versions of SAM: it helps if you delete any existing `SSC`, `SAM`, and `pysam` folders (be sure to commit changes prior to this!). If you have run these scripts previously, also delete the `build_debug`, `build_sam_export`, and `build_ssc_export` folders if they exist. Note that the build scripts below specify a specific tag or branch from each of the projects, modify these as necessary (results not guaranteed...). 

1.  To build SSC and all accompanying projects in their debugging versions, run the following in a command line in the directory `$DEVDIR/neup-ies`:
    ```
    source ./build_debug_SAM
    ```
    - a CodeLite workspace will be created in the directory `$DEVDIR/build_debug`.
    - for more information on debugging with CodeLite, see [these instructions.](https://github.com/uw-esolab/docs/blob/main/sam/debugSSCwithPySSC_Linux_CodeLiteIDE.md) 

2. To build PySAM, everything will be built in their release versions in `build` folders within their respective directories. Run the following in a command line in the directory `$DEVDIR/neup-ies`:
    ```
    source ./build_pysam
    ```
    - a CodeLite workspace will be created in the directories `$DEVDIR/build_ssc_export` and `$DEVDIR/build_sam_export` respectively. 
    - for more information on building PySAM, see [these instructions.](https://github.com/uw-esolab/docs/blob/main/sam/building_PySAM_using_modified_SSC.md)


