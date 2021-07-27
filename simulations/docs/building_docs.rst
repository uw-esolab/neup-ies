.. _building_docs:

Building Documentation
####################################

Documentation here written using Sphinx<http://sphinx-doc.org/>_. 

Make sure you have the proper packages installed in your Python environment::
    
    >>> pip install Sphinx
    >>> pip install numpydoc
    >>> pip install sphinxcontrib-napoleon

While in this directory (``neup-ies/simulations/docs``), make sure you run:

    >>> sphinx-apidoc -M -f -o source/ ../../simulations/
    >>> make html

The final html file will be located at ``neup-ies/simulations/docs/build/html/index.rst``.
