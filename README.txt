
################################################################################
# --- General description

The model simulates spiking networks with different degrees of specific connectivity,
as described in:

Sadeh, Clopath and Rotter (PLOS ONE, 2015).
Processing of Feature Selectivity in Cortical Networks with Specific Connectivity.

Model codes contributed by Sadra Sadeh (s.sadeh@ucl.ac.uk)

Requirements: NEST, Python

[The current codes are written compatible with NEST 2.6.0 and Python 3;
efforts have been made, however, to be backward compatible.]

################################################################################
# --- List of files

[1] SpecNet_source.py
Source class / functions for simulating networks of spiking neurons
with a specified level of specific connectivity in response to oriented stimuli

[2] defaultParams.py
Default parameters for network simulations

[3] SpecNet_run.py
Runs simulations of networks with different degrees of specific connectivity
in response to different stimulus orientations

[4] SpecNet_preprocess
Sample code for preprocessing the raw results of network simulations,
e.g. to extract mean firing rates and tuning curves

################################################################################
# --- Testing the model

(i) Set the parameters of your network simulations in [2];

(ii) Run [3] to simulate the networks and save the resulting simulated data;

(iii) Use [4] to preprocess the raw data and plot example network tuning curves.

################################################################################
