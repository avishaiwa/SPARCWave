####################################
This software package uses both Python and R. All code can and should be run only via Python.

R is used for simulations, wavelet transforms, and the sparcl library (Witten & Tibshirani 2010).
Python is used for CVXPY, and as a wrapper for R code.


Required Python libraries: CVXPY, rpy2, sklearn, joblib, numpy, pylab (and dependencies). 

http://cvxpy.readthedocs.org/en/latest/tutorial/intro/
http://rpy.sourceforge.net/
http://scikit-learn.org/stable/
https://pypi.python.org/pypi/joblib


Required R libraries: sparcl, wmtsa (and dependencies). 
These libraries need to be in the global environment for rpy2 to load them. A simple way to do so is to put
the folders containing these libraries in [R-PATH]/library. For example: "C:\Program Files\R\R-3.2.1\library" 
#####################################

Python files:
SPARCWave_MAIN.py - runs simulations, saves results to disk
SPACWave_functs.py - SPARCWave functions used in other scripts
SPARCWave_get_simulation_results_and_plot.py - loads simulation results from disk, plots
SPARCWave_run_real_data.py - runs SPARCWave on real data files

R files:
R_functions.R - functions for simulations for running the sparcl library

CSV files:
data_growth_T32_Q2_new_ord.csv - Berkeley Growth dataset after scattering transform
labels_growth.csv - labels for Growth dataset

phoneme_scatter_T32_Q8_New_ord.csv - Phoneme dataset after scattering transform
labels_phoneme.csv - labels for Phoneme dataset

wheat_T32_Q8_new_ord.csv - Wheat dataset after scattering transform
labels_wheat.csv - labels for Wheat dataset



