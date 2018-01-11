# Aggregate Size Dependence of Amyloid Adsorption onto Charged Interfaces

This repository contains [Jupyter](http://jupyter.org) Notebooks, as well as experimental and simulation data, for reproducing the work of the scientific paper *Aggregate Size Dependence of Amyloid Adsorption onto Charged Interfaces* published in ACS Langmuir [(doi:10/chz9)](http://pubs.acs.org/doi/abs/10.1021/acs.langmuir.7b03155). This supporting information can be further accessed via Zenodo [(doi:10/ch58)](http://dx.doi.org/10/ch58).

The layout is as follows:

- `notebook.ipynb` - Jupyter Notebook for compiling, running, and plotting all MC simulation data, experimental results, and results from the analytical line segment theory
- `notebook.html` - Jupyter Notebook in HTML format 
- `environment.yml` - Environment file for setting up dependencies for Anaconda
- `mc/` - Directory with Monte Carlo data and C++ code
- `exp/` - Directory with QCM-D and AFM experimental data 
- `md/` - Directory with Notebook for MD simulations using OpenMM and GROMACS input files.

To open this Notebook, install python via [(Mini)conda](https://www.continuum.io/downloads) and make sure all required packages are loaded by issuing the following terminal commands,

    conda env create -f environment.yml
    source activate surfacefibrils
    jupyter-notebook notebook.ipynb
