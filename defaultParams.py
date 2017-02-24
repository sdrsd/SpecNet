################################################################################
# Default parameters for simulating networks of spiking neurons
# with different degrees of specific connectivity
#
# Reference: Sadeh, Clopath and Rotter (PLOS ONE, 2015).
# "Processing of Feature Selectivity in Cortical Networks with Specific Connectivity".
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2015-2016
################################################################################

import numpy as np
import pylab as pl
import os, sys, pickle

#########################
# --- main parameters
#########################

# comment out to have a different randomization of random parameters each time
np.random.seed(1234)

# the range of feature-specific modulations to run the network simulations for
# 0: no modulation in recurrent weights (i.e. unspecific connectivity);
# 1: maximum modulation / specificity of recurrent weights
fs_mod_rng = [0., 0.5, 1.]

# total number of neurons in the network
N = 1000 #5000
# fraction of excitatory neurons
exc_inh = .8
# number of excitatory neurons
NE = int(exc_inh *N)
# number of inhibitory neurons
NI = N - NE

# probability of (E --> E,I) connections
eps_exc = .3
# probability of (I --> E,I) connections
eps_inh = .9

# inhibition dominance of synaptic weights (IPSP = g x EPSP)
g = 10.
# EPSP strength (mV) of recurrent synapses in the local network
J_rec = .3
# (E,I) to (E,I) weights (mV)
J_ee = J_ei = J_rec
J_ie = J_ii = -g*J_rec

# number of neurons in the non-local network
n_ext = 5000.
# spontaneuos firing rate of the non-local network
r_ext = 1.
# overall rate of the non-local input
s_ext = n_ext * r_ext
# EPSP strength (mV) of synapses from non-local sources
J_ext = .2

# number of LGN neurons contacting one cortical neuron
n_ffw = 50.
# average firing rate of LGN neurons
r_ffw = 20.
# overall rate of the feedforward input
r_base = n_ffw * r_ffw
# EPSP strength (mV) of feedforward synapses
J_ffw = 1.

# average synaptic delay (ms)
delay = 1.
# modulation in the input
input_mod = np.array(NE*[0.2] + NI*[.02])

# number of sample (replicate) neurons to record their free membrane potentials
n_smpl = 12

# initial random preferred orientation (PO) of neurons
po_init = np.random.uniform(0., np.pi, N)
po_init[0:n_smpl] = np.arange(0, 180, 180/12.)*np.pi/180
po_init[NE:NE+n_smpl] = np.arange(0, 180, 180/12.)*np.pi/180

# number of stimulus orientations to show
stim_no = 8
# range of stimuli
stim_range = np.arange(0., np.pi, np.pi/stim_no)
# number of trials
trial_no = 1
# simulation time
simtime = 1*1000.
# transient time at the beginning of the simulation
t_trans = 150.

# membrane time constant (ms)
tauMem = 20.
# voltage threshold (mV)
theta = 20.

# neuron parameters (current-based with delta synapses)
neuron_params_delta = {"C_m"      : 1.0,
                "tau_m"     : tauMem,
                "t_ref"     : 2.0,
                "E_L"       : 0.0,
                "V_min"     : -np.inf,
                "V_m"       : 0.,
                "V_reset"   : 0.,
                "V_th"      : theta,
                "I_e"       : 0.}

# parameters of the replicates neuron
# : V_th = inf -> to avoid spiking, to measure the input in absence of spikes
neuron_params_delta_rep = {"C_m"       : 1.0,
                "tau_m"     : tauMem,
                "t_ref"     : 2.0,
                "E_L"       : 0.0,
                "V_min"     : -np.inf,
                "V_m"        : 0.,
                "V_reset"   : 0.,
                "V_th"      : np.inf,
                "I_e"       : 0.}

#########################

# define the NEST path here if needed
nest_path = '/Users/sadra/NEST/nest/ins/lib/python3.4/site-packages/'
if os.path.exists(nest_path): sys.path.append(nest_path)

# code path
code_path = os.getcwd()

# results path
res_path = code_path+'/Results/'
if not os.path.exists(res_path): os.mkdir(res_path)

#########################
#########################
