
################################################################################
# Preprocessing the raw results of network simulations
# to extract firing rates, tuning curves and other desired variables
#
# Reference: Sadeh, Clopath and Rotter (PLOS ONE, 2015).
# "Processing of Feature Selectivity in Cortical Networks with Specific Connectivity".
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2015-2016
################################################################################

import numpy as np; import pylab as pl; import os, time, pickle
import scipy.optimize; from imp import reload
import defaultParams; reload(defaultParams); from defaultParams import *

def sim_results(sim_folder):

    ti = time.time()
    results = {}

    os.chdir(res_path+sim_folder)

    # -- reading spike data and extracting tuning curves
    print('# --- reading spike data\n')

    id_no = n_smpl*4 + 1

    tc_trans, tc = [], []
    for st in enumerate(stim_range):
        print('### stim #: ', str(st[0]+1))

        t0 = time.time()

        tc_st_trans, tc_st = [], []
        for tr in range(trial_no):
            print('# trial #: ', str(tr+1))

            yyy = str(2*N+id_no + tr)
            n0 ='spikes-all-trans-st'+str(st[0])+'-tr'+str(tr)+'-'+yyy+'-0.gdf'
            z = np.loadtxt(n0)

            # tuning curve inferred from transient responses
            fr_trans = np.array([ len(np.where(z[:,0] == n)[0]) / (t_trans/1000) for n in range(1,N+1)])
            tc_st_trans.append(fr_trans)

            yyy = str(2*N+id_no + tr + trial_no)
            n0 ='spikes-all-st'+str(st[0])+'-tr'+str(tr)+'-'+yyy+'-0.gdf'
            z = np.loadtxt(n0)

            # tuning curve exctracted from stationary responses
            fr = np.array([ len(np.where(z[:,0] == n)[0]) / ((simtime-t_trans)/1000) for n in range(1,N+1)])
            tc_st.append(fr)

        t1 = time.time()

        tc.append(tc_st)
        tc_trans.append(tc_st_trans)

    results['tc'] = np.array(tc)
    results['tc_trans'] = np.array(tc_trans)

    # mean tuning curves
    tc_mean = np.mean(tc, 1)

    # saving the results

    f = open('results', 'wb'); pickle.dump(results, f); f.close()

    tf = time.time()

    print('\n# Took : %.2f s' % (tf-ti))

    os.chdir(code_path)

    return tc_mean

# --- Run the preprocessing for all simulations

print('\n##########')
print('Processing data ... \n')

folders = os.listdir(res_path)
TC_all = []
for sim_folder in folders:
    os.chdir(res_path+sim_folder)
    print('Simulation folder --> ', sim_folder, '\n')
    TC = sim_results(sim_folder)
    TC_all.append(TC)
TC_all = np.array(TC_all)

os.chdir(code_path)


# --- Sample figure: plotting network tuning curves
# for different degrees of specific connectivity ---

cls = ['b', 'g', 'r']

pl.figure()

pl.subplot(111); pl.title('Network tuning curves in response to the first orientation')

for i, fs_mod in enumerate(fs_mod_rng):
    pl.plot(po_init[0:NE], TC_all[i,0,0:NE], 'o', mec=cls[i], mfc='none', label=str(fs_mod*100))

pl.legend(title='FS connectivity (%)', loc='best', numpoints=1, frameon=0)

pl.xlabel('Input preferred orientation (radians)', size=15)
pl.ylabel('Firing rate (spikes/s)', size=15)

pl.savefig('SampleNetworkTCs.pdf')

pl.show()
