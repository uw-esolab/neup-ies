# neup-ies

## Integrated Solar & Nuclear Cogeneration of Electricity & Water using the sCO2 Cycle

This is a private repository for sharing code and project files.

*Note from Gabriel: If cloning onto a local machine, I recommend that this repository be cloned into a directory alongside the SAM projects (LK, WEX, SSC, SAM). That is, the **neup-ies** folder should be in the same directory as **lk**, **ssc**, **googletest**, etc. I refer to this directory as `$DEVDIR`.*
<br/><br/>

## Building Documentation

Documentation written using [Sphinx](<http://sphinx-doc.org/>). 
Follow these steps to build the documentation for this project!

Prior to creating the documentation, make sure you have a dedicated Python environment to build everything in. If using [Anaconda](<https://docs.anaconda.com/anaconda/install/linux/>), you could try:

    conda create -n <pysam_env> python=3.7

(without the brackets) or with a name other than `pysam_env`. I have only tested things with Python 3.7, anything else is uncharted territory. Make sure you have the proper packages installed in that Python environment:
    
    pip install Sphinx
    pip install numpydoc
    pip install sphinxcontrib-napoleon

To actually create the docs, make sure you are in the correct conda environment. 
    
    conda activate <pysam_env>
    
Then, run the following:

    cd neup-ies/simulations/docs
    sphinx-apidoc -M -f -o source/ ../../simulations/
    make html

The final html file will be located at ``neup-ies/simulations/docs/build/html/index.html``.

# Running SSC simulations in Linux

Follow the Quick Start Guide in the documentation page to fully build SAM in either Debug or Export (PySAM-compatible) mode. 


