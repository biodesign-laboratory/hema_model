import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
from project_library import lin_sim
from project_library import csv_to_figure_morris
from project_library import merge_figures_grid
import time
import os
import pandas as pd
import M2_debug

# output_names will be used for formatting file names / output_names_laTex will be used when generating matplotlib plots
# same rule applies to param_names / param_names_laTex
# see if-else statements below to determine which preset to use
# this script is specifically meant to be used with MDSC model using morris method

exp_num = 'morris_1'             # used in file names, does not have to be a number

output_names = ['HQ', 'HM', 'N', 'P', 'A', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I']

output_names_laTex = ['H_{Q}', 'H_{M}', 'N', 'P', 'A', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I']

param_names = ['k_H', 'dH', 'theta_N', 'theta_K', 'tau_Q', 'tau_U', 'd_SCSF', 'd_S', 'd_Q', 'd_U', 'd_P', 'd_A', 'g_N', 
                   'N_half', 'S_PH', 'S_PS', 'S_PQ', 'S_AU', 'S_AH', 'S_AS', 'S_AM', 'S_SCSF', 'S_KD', 'k_sn', 'k_nq', 'k_nm', 'k_ns',
                   'R_KU', 'I_crit', 'K_crit', 'psi', 'd_M', 'd_MF', 'S_KMD', 'S_KQ', 'C_QM', 'C_MDM', 'C_UM', 'S_MF', 'omega', 'C_UP', 'alpha']
    
param_names_laTex = ['k_{H}', 'd_{H}', '\\theta_{N}', '\\theta_{K}', '\\tau_{Q}', '\\tau_{U}', 'd_{SCSF}', 'd_{S}', 'd_{Q}', 'd_{U}', 
                    'd_{P}', 'd_{A}', 'g_{N}', 'N_{half}', 'S_{PH}', 'S_{PS}', 'S_{PQ}', 'S_{AU}', 'S_{AH}', 'S_{AS}', 'S_{AM}', 
                    'S_{SCSF}', 'S_{KD}', 'k_{sn}', 'k_{nq}', 'k_{nm}', 'k_{ns}', 'R_{KU}', 'I_{crit}', 'K_{crit}', '\psi', 'd_{M}', 
                    'd_{MF}', 'S_{KMD}', 'S_{KQ}', 'C_{QM}', 'C_{MDM}', 'C_{UM}', 'S_{MF}', '\omega', 'C_{UP}', '\\alpha']

# code directly below checks to see if the relevant folders already exist and if not creates them
script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Model_3_SA')
if not os.path.exists(script_dir):
            os.makedirs(script_dir)

for out_name in output_names:
        
    for order in ['mu', 'mu_star', 'sigma']:

        filepath = os.path.join(script_dir, f'Experiment_{exp_num}', f'{out_name}_out', f'{order}')
        if not os.path.exists(filepath):
            os.makedirs(filepath)

generate_mu_vs_sigma = True
generate_individual_figs = False
generate_merged_figs = False
run_sensitivity_analysis = False
nTimesteps = 100            # how many timesteps to run simulation (default=200)
init_time = 25           # initial time to begin calculating sobol indices, NOT initial time of model(default=100)
delta_t = 1             # timestep size for sobol index calculation, NOT the timestep used in simulation runs

if run_sensitivity_analysis:    # run sensitivity analysis and save relevant sensitivity index to .csv files with delimiter='\t\
    
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
                       [0.007, 0.9],         # d_S
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

    param_values_IHD = morris_sample.sample(problem_IHD, N = 50, num_levels = 8)       # number of samples generated = N * (P + 1) where P = # of parameters   
    print("Shape of the generated sample: ", param_values_IHD.shape)    # sanity check
    print("First few samples:")
    #print(param_values_IHD[:3, ])

    # --- misc code for step size stuff -----
    '''step_sizes = {
    name: (problem_IHD["bounds"][i][1] - problem_IHD["bounds"][i][0]) * (4 / (8 - 1))
    for i, name in enumerate(problem_IHD["names"])
    }
    for name in param_names:
        print(step_sizes[name])'''
    #print(f'Unique values for parameter k_H: {np.unique(param_values_IHD[:, 0])}')
    # -------------------------------

    IHD_out = np.zeros((param_values_IHD.shape[0], 13, int((nTimesteps-init_time)/delta_t)))     # this is where the necessary simulation outputs are stored for calculating indices
    print(f"Shape of array to save outputs to: {IHD_out.shape}")    # sanity check
    # ^ Doesn't store all outputs for memory efficiency purposes

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

    ext_stim = np.zeros((12, int(nTimesteps/0.01)))
    ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

    param_values_dict = [{**dict(zip(param_names, row)), 'N_oo': 2*10**7, 'k': 3} for row in param_values_IHD]           # convert generated samples array into an array of dictionaries

    print("! Computations starting now !")
    start_time = time.time()

    cols = np.arange(int(init_time / 0.01), int(nTimesteps / 0.01), int(delta_t / 0.01))
    print(cols)     # sanity check, make sure this is an array of the desired timesteps to save outputs from

    for i, X in enumerate(param_values_dict):

        output = lin_sim(M2_debug.beta_model_3, X, init_state, nTimesteps, 0.01, ext_stim, ext_stim_m)[0]   # simulate the digital system using sample input

        IHD_out[i] = output[:, cols]
        print(f"Output {i+1} of {param_values_IHD.shape[0]} successfully calculated!")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Outputs generated for all input samples. Elapsed time: ", elapsed_time)

    print("Beginning calculation of morris indices for each timestep now...")    # sanity check
    start_time = time.time()

    for t in range(0, int(nTimesteps - init_time), 1):     # calculate sobol indices from init_time -> end of simulation with 1 timestep in between each

        SI_HQ = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 0, t], num_levels=8)      # args: problem dict, parameter set, model output, num_levels (must be same as num_level in sample method)
        SI_HM = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 1, t], num_levels=8)
        SI_N = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 2, t], num_levels=8)
        SI_P = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 3, t], num_levels=8)
        SI_A = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 4, t], num_levels=8)
        #SI_SCSF = morris_analyze.analyze(problem_IHD, IHD_out[:, 5, t], num_levels=8)
        SI_K = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 6, t], num_levels=8)
        SI_Q = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 7, t], num_levels=8)
        SI_S = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 8, t], num_levels=8)
        SI_U = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 9, t], num_levels=8)
        SI_MDSC = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 10, t], num_levels=8)
        SI_MF = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 11, t], num_levels=8)
        SI_I = morris_analyze.analyze(problem_IHD, param_values_IHD, IHD_out[:, 12, t], num_levels=8)

        #print(SI_HQ.to_df())    # sanity check

        df_SI_HQ = SI_HQ.to_df()
        df_SI_HM = SI_HM.to_df()
        df_SI_N = SI_N.to_df()
        df_SI_P = SI_P.to_df()
        df_SI_A = SI_A.to_df()
        df_SI_K = SI_K.to_df()
        df_SI_Q = SI_Q.to_df()
        df_SI_S = SI_S.to_df()
        df_SI_U = SI_U.to_df()
        df_SI_MDSC = SI_MDSC.to_df()
        df_SI_MF = SI_MF.to_df()
        df_SI_I = SI_I.to_df()

        df_SI_HQ['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HQ_out', 'mu', f'mu_SI_HQ_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_HQ[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HQ_out', 'mu_star', f'mu_star_SI_HQ_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_HQ['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HQ_out', 'sigma', f'sigma_SI_HQ_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        #df_SI_HQ['mu_star_conf'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HQ_out', 'sigma', f'sigma_SI_HQ_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        
        df_SI_HM['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HM_out', 'mu', f'mu_SI_HM_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_HM[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HM_out', 'mu_star', f'mu_star_SI_HM_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_HM['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'HM_out', 'sigma', f'sigma_SI_HM_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_N['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'mu', f'mu_SI_N_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_N[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'mu_star', f'mu_star_SI_N_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_N['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'sigma', f'sigma_SI_N_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_P['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'mu', f'mu_SI_P_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_P[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'mu_star', f'mu_star_SI_P_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_P['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'sigma', f'sigma_SI_P_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_A['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'A_out', 'mu', f'mu_SI_A_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_A[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'A_out', 'mu_star', f'mu_star_SI_A_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_A['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'A_out', 'sigma', f'sigma_SI_A_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        '''mu_SI_SCSF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'SCSF_out', 'mu', f'mu_SI_SCSF_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        mu_star_SI_SCSF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'SCSF_out', 'mu_star', f'mu_star_SI_SCSF_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        sigma_SI_SCSF.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'SCSF_out', 'sigma', f'sigma_SI_SCSF_{exp_num}_{int(t + init_time)}.csv'), sep='\t')'''

        df_SI_K['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'K_out', 'mu', f'mu_SI_K_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_K[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'K_out', 'mu_star', f'mu_star_SI_K_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_K['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'K_out', 'sigma', f'sigma_SI_K_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_S['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'mu', f'mu_SI_S_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_S[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'mu_star', f'mu_star_SI_S_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_S['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'sigma', f'sigma_SI_S_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_Q['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'mu', f'mu_SI_Q_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_Q[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'mu_star', f'mu_star_SI_Q_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_Q['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'sigma', f'sigma_SI_Q_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_U['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'mu', f'mu_SI_U_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_U[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'mu_star', f'mu_star_SI_U_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_U['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'sigma', f'sigma_SI_U_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_MDSC['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MDSC_out', 'mu', f'mu_SI_MDSC_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_MDSC[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MDSC_out', 'mu_star', f'mu_star_SI_MDSC_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_MDSC['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MDSC_out', 'sigma', f'sigma_SI_MDSC_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_MF['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MF_out', 'mu', f'mu_SI_MF_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_MF[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MF_out', 'mu_star', f'mu_star_SI_MF_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_MF['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'MF_out', 'sigma', f'sigma_SI_MF_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        df_SI_I['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'mu', f'mu_SI_I_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_I[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'mu_star', f'mu_star_SI_I_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        df_SI_I['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'sigma', f'sigma_SI_I_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        print(f'Indices for timestep {int(t + init_time)} successfully calculated!')

    end_time = time.time()
    print("All morris indices successfully calculate. Elapsed time: " + str(end_time - start_time))

if generate_individual_figs:    # generate time-series SI graphs for each permutation of (output, parameter)

    print("Figure generation starting now !")
    start_time = time.time()
    orders = ['mu', 'mu_star', 'sigma']

    Filepath = os.path.join(script_dir, f'Experiment_{exp_num}')      # this path should be the branch containing all the output folders which themselves contain the relevant .csv's
    print(f"Filepath: {filepath}")

    for order in orders:

        if order != 'mu_star':
            csv_to_figure_morris(output_names, output_names_laTex, param_names, param_names_laTex, nTimesteps - init_time, init_time, 1, order, exp_num, filepath=Filepath)      # this function takes care of placing the figures in the correct output folders automatically
        else:
            csv_to_figure_morris(output_names, output_names_laTex, param_names, param_names_laTex, nTimesteps - init_time, init_time, 2, order, exp_num, filepath=Filepath, conf_int=True)      # only difference from above is nDatapoints = 2 (because mu_star has conf_int data)
        
        print(f"Figures succesfully generated for {order} indices.")

    
    end_time = time.time()
    print("Individual time-series figures generated for each output and parameter. Elapsed time: " + str(end_time-start_time))

if generate_merged_figs:        # merge time-series SI graphs for each output by parameter into one .png

    output_names = [['HQ', 'HM', 'N', 'P', 'A', 'K'],
                    ['Q', 'S', 'U', 'MDSC', 'MF', 'I']
                    ]

    print("Merged figures generation starting now !")
    start_time = time.time()
    for order in ['mu', 'mu_star', 'sigma']:

        for i in range(2):      # use this to create multiple murged figures, useful when you have many state variables

            merge_figures_grid(2, 3, 800, 600, exp_num, order, output_names[i], param_names, i, script_dir)
            print(f"Merged figures successfully generated for {order} indices.")
    
    end_time = time.time()
    print("Merged figures successfully generated. Elapsed time: " + str(end_time-start_time))

# == Below code is used to plot mu_star against sigma for testing nonlinearity of effects =======
if generate_mu_vs_sigma:

    t = 99          # in this case, we are not plotting multiple timestamps; choose a specific timepoint to evaluate
    Filepath = os.path.join(script_dir, f'Experiment_{exp_num}')


    for i, output in enumerate(output_names):
        
        df_mu_star = pd.read_csv(os.path.join(Filepath, f'{output}_out', 'mu_star', f'mu_star_SI_{output}_{exp_num}_{t}.csv'), delimiter='\t')
        df_sigma = pd.read_csv(os.path.join(Filepath, f'{output}_out', 'sigma', f'sigma_SI_{output}_{exp_num}_{t}.csv'), delimiter='\t')

        df_combined = pd.concat([df_mu_star, df_sigma], axis=1).iloc[:, [0, 1, 4]]

        #print(df_combined.columns)             # sanity check

        plt.scatter(df_combined["mu_star"], df_combined["sigma"], color='blue')
        
        z = max(df_combined['sigma'].max(), df_combined["mu_star"].max())
        if output == 'I':
            print(z)
        plt.plot(np.arange(0, z, 1), np.arange(0, z, 1), 'g--', label = '$f(x)=x$')         # plots f(x) = x
        plt.plot(np.arange(0, z, 1), np.arange(0, 2*z, 2), 'r--', label = '$f(x)=2x$')       # plots f(x) = 2x
        plt.plot(np.arange(0, z, 1), np.arange(0, 1/2*z, 1/2), 'm--', label = '$f(x)=1/2x$')   # plots f(x) = 1/2*x


        for j, row in df_combined.iterrows():
            plt.annotate(row["Unnamed: 0"], (row["mu_star"], row["sigma"]), fontsize=8, xytext=(5,5), textcoords="offset points")

        plt.xlabel("Mu* (Mean Absolute Effect)")
        plt.ylabel("Sigma (Variability)")
        plt.title(f"Mu* vs Sigma (${output_names_laTex[i]}({t}$))")
        plt.legend()
        plt.grid(True)

        path = os.path.join(Filepath, 'mu_star_vs_sigma')
        if not os.path.exists(path):
            os.makedirs(path)

        plt.savefig(os.path.join(path, f'{output}_{exp_num}_{t}.png'))
        plt.close()
        #plt.show()
