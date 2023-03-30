# neup-ies

## Integrated Solar & Nuclear Cogeneration of Electricity & Water using the sCO2 Cycle

This is a repository for sharing code and project files.

*Note from Gabriel: If cloning onto a local machine, I recommend that this repository be cloned into a directory alongside the SAM projects (LK, WEX, SSC, SAM). That is, the **neup-ies** folder should be in the same directory as **lk**, **ssc**, **googletest**, etc. I refer to this directory as `$DEVDIR`.*
<br/><br/>

## Accessing Documentation

You can find a link to the documentation at the following url (hosted by Github Pages):
https://uw-esolab.github.io/neup-ies/

This page is built using .rst files found in the **gh-pages** branch of this project. 
If changes are needed to be made to the documentation, push changes to that branch.

## How to Build Documentation Locally

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

Make sure you are also in the `/neup-ies` directory in your terminal and in the **gh-pages** branch.

Then, run the following:

    sphinx-apidoc -M -f -o docs/_source/ simulations/
    make html

The final html file will be located at ``neup-ies/docs/html/index.html``. This should open in your browser.

## Making Changes to Documentation

If changes made locally to the documentation are ready to be pushed to the main documentation site

    1. switch to the **gh-pages** branch and pull from origin
    2. copy or merge all source rst files to the ``neup-ies/docs/_source/`` directory
    3. copy or merge all built materials to the ``neup-ies/docs/html/`` directory
    4. commit and push changes to remote repository.

Changes may take a couple of minutes to be reflected in the website. 
For reference, I used this guide to help set up the gh-pages: https://python.plainenglish.io/how-to-host-your-sphinx-documentation-on-github-550254f325ae

# Running SSC simulations in Linux

Follow the Quick Start Guide in the documentation page to fully build SAM in either Debug or Export (PySAM-compatible) mode. 


