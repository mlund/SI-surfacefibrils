# Title of manuscript

This repository contain a [Jupyter](http://jupyter.org) Notebook for studying XXX.
The layout is as follows:

- `montecarlo.ipynb` - Jupyter Notebook for compiling, running, and plotting all simulation data
- `environment.yml` - Environment file for setting up dependencies for Anaconda
- `mc/` - Directory with Monte Carlo data and C++ code
- `md/` - Directory with Notebook for MD simulations using OpenMM (experimental!)

To open this Notebook, install python via [(Mini)conda](https://www.continuum.io/downloads) and make sure all required packages are loaded
by issuing the following terminal commands,

    conda env create -f environment.yml
    source activate surfacefibrils
    jupyter-notebook montecarlo.ipynb
