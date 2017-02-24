
################################################################################
# Simulating networks of spiking neurons
# with different degrees of specific connectivity
# in response to different stimulus orientations
#
# Reference: Sadeh, Clopath and Rotter (PLOS ONE, 2015).
# "Processing of Feature Selectivity in Cortical Networks with Specific Connectivity".
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2015-2016
################################################################################

import numpy as np; import pylab as pl; import time, os, pickle
from imp import reload
import defaultParams; reload(defaultParams); from defaultParams import *
import SpecNet_source; reload(SpecNet_source)

# ---
def run_simulation(con_gen, eps_exc, eps_inh, fs_conn, fs_mod):
    '''
    Simulates the response of a network of spiking neurons wit a defined network connectivity
    in response to different stimulus orientations.
    '''

    # instantiation
    MyNet = SpecNet_source.SpecNet()

    # folder
    sim_folder = 'N'+str(N)+ 'g'+str(g)+ '_FS_'+str(fs_conn)+'_Mod'+str(fs_mod)

    if sim_folder in os.listdir(res_path):
       os.chdir(res_path+sim_folder)
    else:
       os.mkdir(res_path+sim_folder)
       os.chdir(res_path+sim_folder)

    # connectivity
    if con_gen:
       # generate the connectivitiy matrix
       con_exc, con_inh = MyNet.random_connectivity(eps_exc, eps_inh)

       con = {}
       con['exc'], con['inh'] = con_exc, con_inh

       # save the connectivitiy matrix
       f = open('con', 'wb'); pickle.dump(con, f); f.close()
    else:
        # load the connectivitiy matrix
        f = open('con', 'rb'); con = pickle.load(f); f.close()

        con_exc, con_inh = con['exc'], con['inh']

    # stimulation
    for stim in enumerate(stim_range):
        print('# ---')
        print('# -- orientation: ', \
        str(stim[0]+1),'/',str(stim_no)+' ('+str(stim[1]*180./np.pi)+' deg)')
        print('# ---')

        stim_id = 'st'+str(stim[0])
        MyNet.stimulate_network(stim_deg=stim[1], stim_id = stim_id, \
        con_exc=con_exc, con_inh=con_inh, mem_pot_rec = 0, fs_conn =fs_conn, fs_mod=fs_mod)

    os.chdir(code_path)

# --- running network simulations for different degrees of specific connectivity
# fs_mod_rng: the range of feature-specific modulations; defined in defaultParams.py
for fs_mod in fs_mod_rng:
    print('\n *** Simulating for feature-specific modulation of: ', fs_mod*100, ' % *** \n')

    # fs_conn: defines whether (E,I) to (E,I) synapses are FS (1) or not (0)
    # order: [EtoE, EtoI, ItoE, ItoI]
    run_simulation(eps_exc=eps_exc, eps_inh=eps_inh, con_gen=1, fs_conn= [1, 0, 0, 0], fs_mod=fs_mod)
