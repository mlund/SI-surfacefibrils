# Aggregate size dependence of amyloid adsorption to charged interfaces

This repository contain a [Jupyter](http://jupyter.org) Notebook for studying amyloid adsorption to an oppositely charged surface.
The layout is as follows:

- `notebook.ipynb` - Jupyter Notebook for compiling, running, and plotting all MD simulation data and results from the analytical line segment theory
- `environment.yml` - Environment file for setting up dependencies for Anaconda
- `mc/` - Directory with Monte Carlo data and C++ code
- `md/` - Directory with Notebook for MD simulations using OpenMM and GROMACS input files.

To open this Notebook, install python via [(Mini)conda](https://www.continuum.io/downloads) and make sure all required packages are loaded
by issuing the following terminal commands,

    conda env create -f environment.yml
    source activate surfacefibrils
    jupyter-notebook notebook.ipynb
