################################################################################
# Source class / functions for simulating networks of spiking neurons
# with different degrees of specific connectivity
#
# Reference: Sadeh, Clopath and Rotter (PLOS ONE, 2015).
# "Processing of Feature Selectivity in Cortical Networks with Specific Connectivity".
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2015-2016
################################################################################

import numpy as np
import pylab as pl
import time, os, sys, pickle
from scipy.stats import norm
from imp import reload
import defaultParams; reload(defaultParams); from defaultParams import *
import nest

################################################################################
## Class for simulating a network of leaky integrate-and-fire (LIF) neurons
################################################################################
class SpecNet(object):

      # --- initial parameters
      def __init__(self):
          pass


      # --- generating random connectivity
      def random_connectivity(self, eps_exc=0.1, eps_inh=0.1):
          '''
          generates random connectivity with defined connection probabilites
          and with fixed in-degrees

          esp_exc: probability of (E --> E,I) connections
          esp_inh: probability of (I --> E,I) connections
          '''

          print("Generating connectivitiy matrix")

          # exc pre-synaptic input per neuron
          CE = int(eps_exc* NE)
          # inh pre-synaptic input per neuron
          CI = int(eps_inh* NI)

          # comment out to generate a new connectivity each time
          np.random.seed(1234)

          # fixed in-degree
          Cin_exc, Cin_inh = [], []
          for n in range(N):
              #if n%1000 == 0: print(n)

              # excitatory
              if n <= NE:
                 zex = list(range(1, n)) + list(range(n+1,NE+1))
              else:
                 zex = list(range(1, NE+1))
              np.random.shuffle(zex)
              cin_exc = zex[0:CE]

              # inhibitory
              if n <= NE:
                 zin = list(range(NE+1, N+1))
              else:
                 zin = list(range(NE+1, n)) + list(range(n+1,N+1))
              np.random.shuffle(zin)
              cin_inh = zin[0:CI]

              Cin_exc.append(cin_exc)
              Cin_inh.append(cin_inh)

          return Cin_exc, Cin_inh


      # --- initiating NEST
      def NEST_initiate(self, dt=.1, n_cores=4):
            '''
            initiating NEST

            dt: time resolution of simulations
            n_cores: how many cores to use for simulations (the more, the faster)
            '''

            # NEST restart
            nest.ResetKernel()
            rng_seed = int(np.random.uniform(0, 12345))
            nest.SetStatus([0],[{'rng_seeds':[rng_seed]}])
            nest.SetKernelStatus({"resolution": dt,
                                "print_time": True,
                                "overwrite_files": True,
                                'local_num_threads': n_cores})
            nest.sr("/RandomConvergentConnect << /allow_multapses false >> SetOptions")


      # --- stimulating the network in response to one stimulus orientation
      def stimulate_network(self, stim_deg, con_exc=[], con_inh=[], stim_id=0, \
                            mem_pot_rec=1, fs_conn = [0,0,0,0], fs_mod=.1):

          # fs_conn: defines whether (E,I) to (E,I) synapses are FS (1) or not (0)
          # order: [EtoE, EtoI, ItoE, ItoI]
          fs_ee, fs_ei, fs_ie, fs_ii = fs_conn

          # --- Start
          ti = time.time()

          self.NEST_initiate()

          print("Network ...")

          # spiking neurons
          nest.SetDefaults("iaf_psc_delta", neuron_params_delta)
          nodes_ex = nest.Create("iaf_psc_delta", NE)
          nodes_in = nest.Create("iaf_psc_delta", NI)

          nodes_all = nodes_ex + nodes_in

          # replicate neurons
          nest.SetDefaults("iaf_psc_delta", neuron_params_delta_rep)
          nodes_ex_rep_excInp = nest.Create("iaf_psc_delta", n_smpl)
          nodes_ex_rep_inhInp = nest.Create("iaf_psc_delta", n_smpl)

          nodes_in_rep_excInp = nest.Create("iaf_psc_delta", n_smpl)
          nodes_in_rep_inhInp = nest.Create("iaf_psc_delta", n_smpl)

          nodes_all_rep_excInp = nodes_ex_rep_excInp + nodes_in_rep_excInp
          nodes_all_rep_inhInp = nodes_ex_rep_inhInp + nodes_in_rep_inhInp

          # vector of total input to all neurons
          input_rate = r_base *(1.+ input_mod* np.cos(2*(stim_deg - po_init)) )

          noise = nest.Create("poisson_generator", N)

          # -- define synapses
          #external
          nest.CopyModel("static_synapse", "ext", {"weight":J_ext, "delay":delay})
          #feedforward
          nest.CopyModel("static_synapse", "ffw", {"weight":J_ffw, "delay":delay})
          # E-to-E
          nest.CopyModel("static_synapse", "ex-ex", {"weight":J_ee, "delay":delay})
          # E-to-I
          nest.CopyModel("static_synapse", "ex-in", {"weight":J_ei, "delay":delay})
          # I-to-E
          nest.CopyModel("static_synapse", "in-ex", {"weight":J_ie, "delay":delay})
          # I-to-I
          nest.CopyModel("static_synapse", "in-in", {"weight":J_ii, "delay":delay})

          # -- devices to record

          # spike detector
          sp_all_trans = nest.Create("spike_detector", trial_no)
          sp_all = nest.Create("spike_detector", trial_no)
          for trial in range(trial_no):
              # transient
              nest.SetStatus([sp_all_trans[trial]],
                          {"label":"spikes-all-trans-"+stim_id+'-tr'+str(trial),
                           "withgid":True, "withtime":True,
                           "to_file":True, "to_memory":False,
                           "start": simtime*trial + 0., "stop": simtime*trial + t_trans })
              nest.ConvergentConnect(nodes_all, [sp_all_trans[trial]], model="ext")
              # stationary
              nest.SetStatus([sp_all[trial]],
                        {"label":"spikes-all-"+stim_id+'-tr'+str(trial),
                         "withgid":True, "withtime":True,
                         "to_file":True, "to_memory":False,
                         "start": simtime*trial + t_trans, "stop": simtime*trial + simtime })
              nest.ConvergentConnect(nodes_all, [sp_all[trial]], model="ext")

          # volt meter
          if mem_pot_rec:
             # stationary
             vm_ex_stat = nest.Create("voltmeter")
             nest.SetStatus(vm_ex_stat, {"label":"vm-exc-"+stim_id,
                         "to_file":True, "to_memory":False})
             nest.DivergentConnect(vm_ex_stat, nodes_ex[0:n_smpl])

             vm_in_stat = nest.Create("voltmeter")
             nest.SetStatus(vm_in_stat, {"label":"vm-inh-"+stim_id,
                         "to_file":True, "to_memory":False})
             nest.DivergentConnect(vm_in_stat, nodes_in[0:n_smpl])

             # replicates (measuring input)
             vm_ex_excInp = nest.Create("voltmeter")
             nest.SetStatus(vm_ex_excInp, {"label":"vm-exc-excInp-"+stim_id,
                         "to_file":True, "to_memory":False})
             nest.DivergentConnect(vm_ex_excInp, nodes_ex_rep_excInp[0:n_smpl])

             vm_ex_inhInp = nest.Create("voltmeter")
             nest.SetStatus(vm_ex_inhInp, {"label":"vm-exc-inhInp-"+stim_id,
                         "to_file":True, "to_memory":False})
             nest.DivergentConnect(vm_ex_inhInp, nodes_ex_rep_inhInp[0:n_smpl])

             vm_in_excInp = nest.Create("voltmeter")
             nest.SetStatus(vm_in_excInp, {"label":"vm-inh-excInp-"+stim_id,
                         "to_file":True, "to_memory":False})
             nest.DivergentConnect(vm_in_excInp, nodes_in_rep_excInp[0:n_smpl])

             vm_in_inhInp = nest.Create("voltmeter")
             nest.SetStatus(vm_in_inhInp, {"label":"vm-inh-inhInp-"+stim_id,
                         "to_file":True, "to_memory":False})
             nest.DivergentConnect(vm_in_inhInp, nodes_in_rep_inhInp[0:n_smpl])

          # -- Connecting the network
          print("Connecting the network ...")

          CE = int(eps_exc* NE)
          CI = int(eps_inh* NI)

          for nn in enumerate(nodes_all):
              nest.SetStatus([noise[nn[0]]], {'rate':input_rate[nn[0]]})
              nest.Connect([noise[nn[0]]], [nn[1]], model="ffw")

          ext_inp = nest.Create("poisson_generator")
          nest.SetStatus(ext_inp, {"rate":r_ext})
          nest.DivergentConnect(ext_inp, nodes_all, model="ext")

          # -- feature-specific modulation of the weights
          for ne in range(NE):
              # exc --> exc
               dth = po_init[np.array(con_exc[ne])-1] - po_init[ne]
               exc_list = (J_ee *(1.+ fs_ee *fs_mod * np.cos(2*dth))).tolist()
               nest.ConvergentConnect(con_exc[ne], [nodes_all[ne]], exc_list, CE*[delay])
               if ne < n_smpl:
                   nest.ConvergentConnect(con_exc[ne], [nodes_all_rep_excInp[ne]], exc_list, CE*[delay])
              # inh --> exc
               dth = po_init[np.array(con_inh[ne])-1] - po_init[ne]
               inh_list = (J_ie *(1+fs_ei*fs_mod *np.cos(2*dth))).tolist()
               nest.ConvergentConnect(con_inh[ne], [nodes_all[ne]], inh_list, CI*[delay])
               if ne < n_smpl:
                   nest.ConvergentConnect(con_inh[ne], [nodes_all_rep_inhInp[ne]], inh_list, CI*[delay])

          for ni in range(NE, N):
              # exc --> inh
               dth = po_init[np.array(con_exc[ni])-1] - po_init[ni]
               exc_list = (J_ei *(1+fs_ei*fs_mod *np.cos(2*dth))).tolist()
               nest.ConvergentConnect(con_exc[ni], [nodes_all[ni]], exc_list, CE*[delay])
              # inh --> inh
               dth = po_init[np.array(con_inh[ni])-1] - po_init[ni]
               inh_list = (J_ii *(1+fs_ii*fs_mod *np.cos(2*dth))).tolist()
               nest.ConvergentConnect(con_inh[ni], [nodes_all[ni]], inh_list, CI*[delay])

          # --- running simulations

          print("Simulating the network ...")

          for trial in range(trial_no):
              print('# -- Trial # ', str(trial+1))
              nest.Simulate(simtime + t_trans)

          ts = time.time()
          sim_time = ts - ti
          r_avg = nest.GetStatus([sp_all[0]], 'n_events')[0]/ (N*simtime/1000.)

          print('\n########################################')
          print('End of simulation.')
          print("Simulation time   : %.2f s" % sim_time)
          print("Average rate      : %.2f Hz" %  r_avg)
          print('######################################## \n')
