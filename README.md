# Coupled-Cluster-NuclearMatter

The coupled-cluster calculations (f90 files) can be compiled using mpifort linked with LAPACK libraries and openMP and run
using:  mpirun -np <#processes> ccm.exe ccm_in
the ccm_in file is the input file that specifies the momentum grid, number of particles, solution method (e.g., include triples), etc...

The data-analysis folder contains a python script to plot some example data using a bayesian formula to approx eft error bars.
The NeuralNetworkRegression.py file was a preliminary attempt to train a artificial NN to predict calculations given several input parameter to ultimately aid in fitting.
