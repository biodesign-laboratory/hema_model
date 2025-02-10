import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from SALib.test_functions import Ishigami
from project_library import lin_sim
from project_library import csv_to_figure
from project_library import merge_figures_grid
import time
import pandas
import os
import M2_debug

# output_names will be used for formatting file names / output_names_laTex will be used when generating matplotlib plots
# same rule applies to param_names / param_names_laTex
# see if-else statements below to determine which preset to use

# ========= commented out below, use this script only for model 3 =====================
'''output_preset = 2
param_preset = 1

if output_preset == 1:      # H(t), S(t), I(t), Q(t), U(t)

    output_names = ['H', 'S', 'I', 'N', 'Q', 'U']
    output_names_laTex = ['H', 'S', 'I', 'N', 'Q', 'U']

elif output_preset == 2:    # HQ(t), HM(t), N(t), P(t), A(t), SCSF(t), K(t), Q(t), S(t), U(t), MDSC(t), MF(t), I(t); use for model 3 (MDSC model)

    output_names = ['HQ', 'HM', 'N', 'P', 'A', 'SCSF', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I']
    output_names_laTex = ['H_{Q}', 'H_{M}', 'N', 'P', 'A', 'SCSF', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I']


if param_preset == 1:      # test all derivate parameters as well as H_init, N_init

    param_names = ['gn', 'knq', 'kns', 'ktn', 'w', 'pcrit', 'spq', 'sph',
                'sps', 'sas', 'sah', 'sau', 'Is', 'As', 'Ic', 'Ac',
                'y', 'ds', 'dp', 'da', 'dq', 'du', 'H_init', 'N_init']
    
    param_names_laTex = ['g_{n}', 'k_{nq}', 'k_{ns}', 'k_{tn}', 'w', 'p_{crit}', 's_{pq}', 's_{ph}',
                's_{ps}', 's_{as}', 's_{ah}', 's_{au}', 'I_{s}', 'A_{s}', 'I_{c}', 'A_{c}',
                'y', 'd_{s}', 'd_{p}', 'd_{a}', 'd_{q}', 'd_{u}', 'H_{Init}', 'N_{Init}']
    
elif param_preset == 2:    # same as 1 but doesn't contain N_init

    param_names = ['gn', 'knq', 'kns', 'ktn', 'w', 'pcrit', 'spq', 'sph',
                'sps', 'sas', 'sah', 'sau', 'Is', 'As', 'Ic', 'Ac',
                'y', 'ds', 'dp', 'da', 'dq', 'du', 'H_init']

    param_names_laTex = ['g_{n}', 'k_{nq}', 'k_{ns}', 'k_{tn}', '\omega', 'p_{crit}', 's_{pq}', 's_{ph}',
                's_{ps}', 's_{as}', 's_{ah}', 's_{au}', 'I_{s}', 'A_{s}', 'I_{c}', 'A_{c}',
                '\gamma', 'd_{s}', 'd_{p}', 'd_{a}', 'd_{q}', 'd_{u}', 'H_{Init}']
    
elif param_preset == 3:     # test all parameters for model 3

    param_names = ['k_H', 'dH', 'theta_N', 'theta_K', 'tau_Q', 'tau_U', 'dSCSF', 'd_S', 'd_Q', 'd_U', 'd_P', 'd_A', 'g_N', 
                   'N_half', 'S_PH', 'S_PS', 'S_PQ', 'S_AU', 'S_AH', 'S_AS', 'S_AM', 'S_SCSF', 'S_KD', 'k_sn', 'k_nq', 'k_nm', 'k_ns',
                   'R_KU', 'I_crit', 'K_crit', 'k', 'psi', 'd_M', 'S_KMD', 'S_KQ', 'C_QM', 'C_MDM', 'C_UM', 'S_MF']
    
    param_names_laTex = ['k_{H}', 'd_{H}', '\theta_{N}', '\theta_{K}', '\tau_{Q}', '\tau_{U}', 'd_{SCSF}', 'd_{S}', 'd_{Q}', 'd_{U}', 
                         'd_{P}', 'd_{A}', 'g_{N}', 'N_{half}', 'S_{PH}', 'S_{PS}', 'S_{PQ}', 'S_{AU}', 'S_{AH}', 'S_{AS}', 'S_{AM}', 
                         'S_{SCSF}', 'S_{KD}', 'k_{sn}', 'k_{nq}', 'k_{nm}', 'k_{ns}', 'R_{KU}', 'I_{crit}', 'K_{crit}', 'k', '\psi', 'd_{M}', 
                         'd_{MF}', 'S_{KMD}', 'S_{KQ}', 'C_{QM}', 'C_{MDM}', 'C_{UM}', 'S_{MF}']'''
    
# ========= commented out above, use this script only for model 3 =====================

output_names = ['HQ', 'HM', 'N', 'P', 'A', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I']

output_names_laTex = ['H_{Q}', 'H_{M}', 'N', 'P', 'A', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I']

param_names = ['k_H', 'dH', 'theta_N', 'theta_K', 'tau_Q', 'tau_U', 'd_SCSF', 'd_S', 'd_Q', 'd_U', 'd_P', 'd_A', 'g_N', 
                   'N_half', 'S_PH', 'S_PS', 'S_PQ', 'S_AU', 'S_AH', 'S_AS', 'S_AM', 'S_SCSF', 'S_KD', 'k_sn', 'k_nq', 'k_nm', 'k_ns',
                   'R_KU', 'I_crit', 'K_crit', 'psi', 'd_M', 'd_MF', 'S_KMD', 'S_KQ', 'C_QM', 'C_MDM', 'C_UM', 'S_MF', 'omega', 'C_UP', 'alpha']
    
param_names_laTex = ['k_{H}', 'd_{H}', '\\theta_{N}', '\\theta_{K}', '\\tau_{Q}', '\\tau_{U}', 'd_{SCSF}', 'd_{S}', 'd_{Q}', 'd_{U}', 
                    'd_{P}', 'd_{A}', 'g_{N}', 'N_{half}', 'S_{PH}', 'S_{PS}', 'S_{PQ}', 'S_{AU}', 'S_{AH}', 'S_{AS}', 'S_{AM}', 
                    'S_{SCSF}', 'S_{KD}', 'k_{sn}', 'k_{nq}', 'k_{nm}', 'k_{ns}', 'R_{KU}', 'I_{crit}', 'K_{crit}', '\psi', 'd_{M}', 
                    'd_{MF}', 'S_{KMD}', 'S_{KQ}', 'C_{QM}', 'C_{MDM}', 'C_{UM}', 'S_{MF}', '\omega', 'C_{UP}', '\\alpha']

# What is absent from parameters being generated?
#   - N_oo (fixed at 2 * 10**7 (20 million))
#   - k (hill-type coefficient, fixed at 3, mostly arbitrary)

exp_num = 1             # trial number
generate_individual_figs = True
generate_merged_figs = True
run_sensitivity_analysis = True
nTimesteps = 100            # how many timesteps to run simulation (default=200)
init_time = 75           # initial time to begin calculating sobol indices, NOT initial time of model(default=100)
#t_final = 200
#model_to_test = 2       # 1 -> no MDSCs (Model 2), 2 -> MDSCs (Model 3); relevant only if run_sensitivity_analysis = True; ignore this

if run_sensitivity_analysis:    # run sensitivity analysis and save relevant sensitivity index to .csv files with delimiter='\t\

    '''if param_preset == 1:

        problem_IHD = {
            'num_vars': 24,
            'names': ['g_n', 'k_nq', 'k_ns', 'k_tn', 'w', 'p_crit', 's_pq', 's_ph', 's_ps', 's_as', 's_ah', 's_au', 
                    'I_s', 'A_s', 'I_c', 'A_c', 'y', 'd_s', 'd_p', 'd_a', 'd_q', 'd_u', 'H_init', 'N_init'],
            'bounds': [[0.5, 2.5],
                    [3, 10],
                    [1, 3],
                    [0.5, 2],
                    [0.00005, 0.005],
                    [3000, 10000],
                    [2, 5],
                    [0.25, 1],
                    [0.5, 2],
                    [0.5, 2],
                    [0.25, 1],
                    [2, 5],
                    [1000, 5000],
                    [2000, 6000],
                    [4000, 10000],
                    [3000, 8000],
                    [0.00005, 0.005],
                    [0.05, 0.25],
                    [0.75, 1.5],
                    [0.5,1.25],
                    [0.75, 2],
                    [0.5, 1.5],
                    [500, 5000],
                    [0, 20000]]
        }

    elif param_preset == 2:
        
        problem_IHD = {
            'num_vars': 23,
            'names': ['g_n', 'k_nq', 'k_ns', 'k_tn', 'w', 'p_crit', 's_pq', 's_ph', 's_ps', 's_as', 's_ah', 's_au', 
                    'I_s', 'A_s', 'I_c', 'A_c', 'y', 'd_s', 'd_p', 'd_a', 'd_q', 'd_u', 'H_init'],
            'bounds': [[0.5, 2.5],
                    [3, 10],
                    [1, 3],
                    [0.5, 2],
                    [0.00005, 0.005],
                    [3000, 10000],
                    [2, 5],
                    [0.25, 1],
                    [0.5, 2],
                    [0.5, 2],
                    [0.25, 1],
                    [2, 5],
                    [1000, 5000],
                    [2000, 6000],
                    [4000, 10000],
                    [3000, 8000],
                    [0.00005, 0.005],
                    [0.05, 0.25],
                    [0.75, 1.5],
                    [0.5,1.25],
                    [0.75, 2],
                    [0.5, 1.5],
                    [500, 20000]]
        }
    
    elif param_preset == 3:
        
        problem_IHD = {
            'num_vars': 29,
            'names': ['g_n', 'k_nq', 'k_ns', 'k_tn', 'w', 'p_crit', 's_pq', 's_ph', 's_ps', 's_as', 's_ah', 's_au', 
                    'I_s', 'A_s', 'I_c', 'A_c', 'y', 'd_s', 'd_p', 'd_a', 'd_q', 'd_u', 'H(t)', 'N(t)', 'P(t)', 'A(t)', 'Q(t)', 'S(t)', 'U(t)', 'H_init'],
            'bounds': [[0.5, 2.5],              # g_n
                    [3, 20],                    # k_nq
                    [1, 5],                     # k_ns
                    [0.5, 10],                  # k_tn
                    [0.000005, 0.05],           # w
                    [1000, 15000],              # p_crit
                    [2, 10],                    # s_pq
                    [0.25, 3],                  # s_ph
                    [0.5, 5],                   # s_ps
                    [0.5, 5],                   # s_as
                    [0.25, 3],                  # s_ah
                    [2, 10],                    # s_au
                    [1000, 10000],              # I_s
                    [1000, 10000],              # A_s
                    [4000, 20000],              # I_c
                    [4000, 20000],              # A_c
                    [0.000005, 0.05],           # y
                    [0.05, 0.5],                # d_s
                    [0.2, 1.5],                # d_p
                    [0.2, 1.5],                # d_a
                    [0.2, 1.5],                # d_q
                    [0.2, 1.5],                 # d_u
                    [0, 20000],                 # H(t)
                    [0, 20000],                 # N(t)
                    [0, 20000],                 # P(t)
                    [0, 20000],                 # A(t) 
                    [0, 20000],                 # Q(t)
                    [0, 20000],                 # S(t)
                    [0, 20000],                 # U(t)
                    [0, 20000]]                 # H_init
        }'''
    
    problem_IHD = {
            'num_vars': 42,
            'names': ['k_H', 'dH', 'theta_N', 'theta_K', 'tau_Q', 'tau_U', 'dSCSF', 'd_S', 'd_Q', 'd_U', 'd_P', 'd_A', 'g_N', 
                   'N_half', 'S_PH', 'S_PS', 'S_PQ', 'S_AU', 'S_AH', 'S_AS', 'S_AM', 'S_SCSF', 'S_KD', 'k_sn', 'k_nq', 'k_nm', 'k_ns',
                   'R_KU', 'I_crit', 'K_crit', 'psi', 'd_M', 'd_MF', 'S_KMD', 'S_KQ', 'C_QM', 'C_MDM', 'C_UM', 'S_MF', 'omega', 'C_UP', 'alpha'],
            'bounds': [[1, 10],         # k_H
                       [0.05, 0.5],         # dH
                       [5000, 50_000],         # theta_N
                       [5000, 50_000],         # theta_K
                       [1, 10],         # tau_Q
                       [1, 10],         # tau_U
                       [0.05, 0.5],         # d_SCSF
                       [0.007, 0.1],         # d_S
                       [0.7, 1],         # d_Q
                       [0.25, 0.85],         # d_U
                       [0.7, 1],         # d_P
                       [0.7, 1],         # d_A
                       [0.002, 0.2],         # g_N
                       [50, 5000],         # N_half
                       [1, 10],         # S_PH
                       [1, 10],         # S_PS
                       [5, 50],         # S_PQ
                       [5, 50],         # S_AU
                       [1, 10],         # S_AH
                       [1, 10],         # S_AS
                       [3, 30],         # S_AM
                       [1000, 20000],         # S_SCSF
                       [1, 10],         # S_KD
                       [1, 10],         # k_sn
                       [2, 20],         # k_nq
                       [1, 10],         # k_nm
                       [1, 10],         # k_ns
                       [1, 100],         # R_KU
                       [0.6, 0.95],         # I_crit
                       [2000, 40_000],         # K_crit
                       [0.1, 1],         # psi
                       [0.5, 1],         # d_M
                       [0.2, 0.9],         # d_MF
                       [3, 30],         # S_KMD
                       [5, 50],         # S_KQ
                       [1, 10],         # C_QM
                       [1, 10],         # C_MDM
                       [1, 10],         # C_UM
                       [1000, 10_000],         # S_MF
                       [0.05, 0.95],    # omega
                       [1, 10],         # C_UP
                       [0.1, 1]         # alpha
                       ]                
        }

    param_values_IHD = sobol_sample.sample(problem_IHD, 2048, calc_second_order=True)      # use n=2 when checking for errors, else use n=2048
    print("Shape of the generated sample: ", param_values_IHD.shape)
    print("First few samples:")
    print(param_values_IHD[:3])

    '''init_values = np.hstack((param_values_IHD[:, -2:-1], np.tile(np.array([0, 4000, 0, 0.9, 0.1, 1800, 300, 1]), (param_values_IHD.shape[0], 1))))
    print("First few initial value sets: ")         # sanity check
    print(init_values[3,:])
    print("Shape of init_value array: ", init_values.shape)'''

    IHD_out = np.zeros((param_values_IHD.shape[0], 13, int(nTimesteps/0.01)))     # this is where the simulation outputs are stored

    init_state = [

        10000,  # Quiescent HSPCs
        5000,      # Proliferating HSPCs
        0,      # PAMPs (Pathogens)
        10000,      # Pro-inflammatory Cytokines
        10000,      # Anti-inflammatory Cytokines
        10000,  # Stem Cell Supporting Factors
        0,      # DAMPs (Tissue Damage)
        500,      # Activated leukocytes
        500,   # Stable leukocytes
        500,       # Suppressor leukocytes
        1,      # MDSC
        5000,  # MF

    ]       # init_state is mostly arbitary, need to check how this affects interpretability of SA results

    ext_stim = np.zeros((param_values_IHD.shape[0], 12, int(nTimesteps/0.01)))
    ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

    param_values_dict = [{**dict(zip(param_names, row)), 'N_oo': 2*10**7, 'k': 3} for row in param_values_IHD]           # convert generated samples array into an array of dicts

    print("! Computations starting now !")
    start_time = time.time()
    for i, X in enumerate(param_values_dict):

        output = lin_sim(M2_debug.beta_model_3, X, init_state, nTimesteps, 0.01, ext_stim[i], ext_stim_m)[0]   # simulate the digital system using sample input

        '''for x in range(init_time, init_time + 100):

            # IHD_out[i, :, x] = output[:, int(x/0.01)]'''
        '''print(IHD_out[i].shape)
        print(output.shape)'''
        IHD_out[i] = output

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Outputs generated for all input samples. Elapsed time: ", elapsed_time)

    print("Beginning calculation of sobol indices for each timestep now...")    # sanity check
    start_time = time.time()

    # code directly below checks to see if the relevant folders already exist and if not creates them
    script_dir = os.path.dirname(os.path.abspath(__file__))

    for out_name in output_names:
            
        for order in ['first', 'second', 'total']:

            filepath = os.path.join(script_dir, 'Model_3_SA', f'Experiment_{exp_num}', f'{out_name}_out', f'{order}_order')
            if not os.path.exists(filepath):
                os.makedirs(filepath)


    for t in range(int(init_time/0.01), int(nTimesteps/0.01), 100):     # calculate sobol indices from init_time -> end of simulation with 1 timestep in between each
        
        # print statements below are sanity checks
        '''print("Param values shape: ", param_values_IHD.shape)
        #print(IHD_out[:, 0, t].shape)

        if np.any(np.isnan(IHD_out[:, 0, t])) or np.any(np.isinf(IHD_out[:, 0, t])):
            raise ValueError("IHD_out contains NaN or Inf values!")
        
        #output_fake = np.random.rand(param_values_IHD.shape[0])        # for debugging purposes

        #print("problem_IHD:", problem_IHD)
        print("Number of variables:", problem_IHD['num_vars'])
        print("Bounds length:", len(problem_IHD['bounds']))
        #print("Bounds:", problem_IHD['bounds'])
        print("Param values [0] shape: ", param_values_IHD.shape[0])
        print("Output[0] shape: ", IHD_out[:, 0, t].shape)'''

        SI_HQ = sobol_analyze.analyze(problem_IHD, IHD_out[:, 0, t], calc_second_order=True)
        SI_HM = sobol_analyze.analyze(problem_IHD, IHD_out[:, 1, t], calc_second_order=True)
        SI_N = sobol_analyze.analyze(problem_IHD, IHD_out[:, 2, t], calc_second_order=True)
        SI_P = sobol_analyze.analyze(problem_IHD, IHD_out[:, 3, t], calc_second_order=True)
        SI_A = sobol_analyze.analyze(problem_IHD, IHD_out[:, 4, t], calc_second_order=True)
        #SI_SCSF = sobol_analyze.analyze(problem_IHD, IHD_out[:, 5, t], calc_second_order=True)
        SI_K = sobol_analyze.analyze(problem_IHD, IHD_out[:, 6, t], calc_second_order=True)
        SI_Q = sobol_analyze.analyze(problem_IHD, IHD_out[:, 7, t], calc_second_order=True)
        SI_S = sobol_analyze.analyze(problem_IHD, IHD_out[:, 8, t], calc_second_order=True)
        SI_U = sobol_analyze.analyze(problem_IHD, IHD_out[:, 9, t], calc_second_order=True)
        SI_MDSC = sobol_analyze.analyze(problem_IHD, IHD_out[:, 10, t], calc_second_order=True)
        SI_MF = sobol_analyze.analyze(problem_IHD, IHD_out[:, 11, t], calc_second_order=True)
        SI_I = sobol_analyze.analyze(problem_IHD, IHD_out[:, 12, t], calc_second_order=True)

        total_SI_HQ, first_SI_HQ, second_SI_HQ = SI_HQ.to_df()
        total_SI_HM, first_SI_HM, second_SI_HM = SI_HM.to_df()
        total_Si_N, first_Si_N, second_Si_N = SI_N.to_df()
        total_Si_P, first_Si_P, second_Si_P = SI_P.to_df()
        total_Si_A, first_Si_A, second_Si_A = SI_A.to_df()
        #total_Si_SCSF, first_Si_SCSF, second_Si_SCSF = SI_SCSF.to_df()
        total_Si_K, first_Si_K, second_Si_K = SI_K.to_df()
        total_Si_S, first_Si_S, second_Si_S = SI_S.to_df()
        total_Si_Q, first_Si_Q, second_Si_Q = SI_Q.to_df()
        total_Si_U, first_Si_U, second_Si_U = SI_U.to_df()
        total_Si_MDSC, first_Si_MDSC, second_Si_MDSC = SI_MDSC.to_df()
        total_Si_MF, first_Si_MF, second_Si_MF = SI_U.to_df()
        total_Si_I, first_Si_I, second_Si_I = SI_I.to_df()

        total_SI_HQ.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HQ_out', 'total_order', f'total_SI_HQ_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_SI_HQ.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HQ_out', 'first_order', f'first_SI_HQ_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_SI_HQ.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HQ_out', 'second_order', f'second_SI_HQ_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        
        total_SI_HM.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HM_out', 'total_order', f'total_SI_HM_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_SI_HM.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HM_out', 'first_order', f'first_SI_HM_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_SI_HM.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HM_out', 'second_order', f'second_SI_HM_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'total_order', f'total_Si_N_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'first_order', f'first_Si_N_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'second_order', f'second_Si_N_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_P.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'total_order', f'total_Si_P_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_P.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'first_order', f'first_Si_P_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_P.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'second_order', f'second_Si_P_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_A.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'A_out', 'total_order', f'total_Si_A_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_A.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'A_out', 'first_order', f'first_Si_A_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_A.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'A_out', 'second_order', f'second_Si_A_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        '''total_Si_SCSF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'SCSF_out', 'total_order', f'total_Si_SCSF_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_SCSF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'SCSF_out', 'first_order', f'first_Si_SCSF_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_SCSF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'SCSF_out', 'second_order', f'second_Si_SCSF_{exp_num}_{int(t*0.01)}.csv'), sep='\t')'''

        total_Si_K.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'K_out', 'total_order', f'total_Si_K_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_K.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'K_out', 'first_order', f'first_Si_K_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_K.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'K_out', 'second_order', f'second_Si_K_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_S.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'total_order', f'total_Si_S_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_S.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'first_order', f'first_Si_S_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_S.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'second_order', f'second_Si_S_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_Q.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'total_order', f'total_Si_Q_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_Q.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'first_order', f'first_Si_Q_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_Q.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'second_order', f'second_Si_Q_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_U.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'total_order', f'total_Si_U_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_U.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'first_order', f'first_Si_U_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_U.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'second_order', f'second_Si_U_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_MDSC.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MDSC_out', 'total_order', f'total_Si_MDSC_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_MDSC.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MDSC_out', 'first_order', f'first_Si_MDSC_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_MDSC.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MDSC_out', 'second_order', f'second_Si_MDSC_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_MF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MF_out', 'total_order', f'total_Si_MF_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_MF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MF_out', 'first_order', f'first_Si_MF_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_MF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MF_out', 'second_order', f'second_Si_MF_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        total_Si_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'total_order', f'total_Si_I_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        first_Si_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'first_order', f'first_Si_I_{exp_num}_{int(t*0.01)}.csv'), sep='\t')
        second_Si_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'second_order', f'second_Si_I_{exp_num}_{int(t*0.01)}.csv'), sep='\t')

        print(f'Indices for timestep {int(t*0.01)} successfully calculated!')

    end_time = time.time()
    print("All sobol indices successfully calculate. Elapsed time: " + str(end_time - start_time))

if generate_individual_figs:    # generate time-series SI graphs for each permutation of (output, parameter)

    print("Figure generation starting now !")
    start_time = time.time()
    orders = ['first', 'total']

    for order in orders:

        filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'Experiment_{exp_num}')      # this path should be the branch containing all the output folders which themselves contain the relevant .csv's
        csv_to_figure(output_names, output_names_laTex, param_names, param_names_laTex, nTimesteps - init_time, init_time, 2, order, exp_num, filepath)      # this function takes care of placing the figures in the correct output folders automatically

        print(f"Figures succesfully generated for {order}-order indices.")

    
    end_time = time.time()
    print("Individual time-series figures generated for each output and parameter. Elapsed time: " + str(end_time-start_time))

if generate_merged_figs:        # merge time-series SI graphs for each output by parameter into one .png

    output_names = [['HQ', 'HM', 'N', 'P', 'A', 'K'],
                    ['Q', 'S', 'U', 'MDSC', 'MF', 'I']
                    ]

    print("Merged figures generation starting now !")
    start_time = time.time()
    for order in ['first', 'total']:

        for i in range(2):      # use this to create multiple murged figures, useful when you have many state variables

            merge_figures_grid(2, 3, 800, 600, exp_num, order, output_names[i], param_names, i)
            print(f"Merged figures successfully generated for {order}-order indices.")

        '''else:
            merge_figures_grid(2, 3, 800, 600, exp_num, order, output_names, param_names)
            print(f"Merged figures successfully generated for {order}-order indices.")'''           # else don't change output_names set at beginning and uncomment this
    
    end_time = time.time()
    print("Merged figures successfully generated. Elapsed time: " + str(end_time-start_time))\
    
# note that second_order figures must be handled separately due to the difference in data shape
