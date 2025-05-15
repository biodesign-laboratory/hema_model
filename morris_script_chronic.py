import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
from project_library import lin_sim
from project_library import chronic_infection_flagger
import time
import os
import pandas as pd
import M3_beta

# output_names will be used for formatting file names / output_names_laTex will be used when generating matplotlib plots
# same rule applies to param_names / param_names_laTex
# see if-else statements below to determine which preset to use

exp_num = 'name_here'             # used in file names, does not have to be a number

output_names = ['chronic_duration']

output_names_laTex = ['Chronic Duration']

param_names = ['k_H', 'dH', 'theta_N', 'theta_K', 'tau_Q', 'tau_U', 'd_SCSF', 'd_S', 'd_Q', 'd_U', 
                   'N_half', 'S_PH', 'S_PS', 'S_PQ', 'S_AU', 'S_AH', 'S_AS', 'S_AM', 'S_SCSF', 'S_KD', 'k_sn', 'k_nq', 'k_nm', 'k_ns',
                   'R_KU', 'I_crit', 'K_crit', 'psi', 'd_M', 'd_MF', 'S_KMD', 'S_KQ', 'C_QM', 'C_MDM', 'C_UM', 'S_MF', 'omega', 'C_UP', 'alpha', 'H_crit', 'pathogen_input']
    
param_names_laTex = ['k_{H}', 'd_{H}', '\\theta_{N}', '\\theta_{K}', '\\tau_{Q}', '\\tau_{U}', 'd_{SCSF}', 'd_{S}', 'd_{Q}', 'd_{U}',
                    'N_{half}', 'S_{PH}', 'S_{PS}', 'S_{PQ}', 'S_{AU}', 'S_{AH}', 'S_{AS}', 'S_{AM}', 
                    'S_{SCSF}', 'S_{KD}', 'k_{sn}', 'k_{nq}', 'k_{nm}', 'k_{ns}', 'R_{KU}', 'I_{crit}', 'K_{crit}', '\psi', 'd_{M}', 
                    'd_{MF}', 'S_{KMD}', 'S_{KQ}', 'C_{QM}', 'C_{MDM}', 'C_{UM}', 'S_{MF}', '\omega', 'C_{UP}', '\\alpha', 'H_{crit}', 'P_{in}']

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
run_sensitivity_analysis = True
nTimesteps = 400            # how many timesteps to run simulation (default=200)
# init_time = 25           # initial time to begin calculating sobol indices, NOT initial time of model(default=100)

if run_sensitivity_analysis:    # run sensitivity analysis and save relevant sensitivity index to .csv files with delimiter='\t\
    
    problem_IHD = {
            'num_vars': 41,
            'names': ['k_H', 'dH', 'theta_N', 'theta_K', 'tau_Q', 'tau_U', 'dSCSF', 'd_S', 'd_Q', 'd_U', 
                   'N_half', 'S_PH', 'S_PS', 'S_PQ', 'S_AU', 'S_AH', 'S_AS', 'S_AM', 'S_SCSF', 'S_KD', 'k_sn', 'k_nq', 'k_nm', 'k_ns',
                   'R_KU', 'I_crit', 'K_crit', 'psi', 'd_M', 'd_MF', 'S_KMD', 'S_KQ', 'C_QM', 'C_MDM', 'C_UM', 'S_MF', 'omega', 'C_UP', 'alpha', 'H_crit', 'path_size'],
            'bounds': [[1, 10],                 # k_H
                       [0.05, 0.8],             # dH
                       [5000, 50_000],          # theta_N
                       [5000, 50_000],          # theta_K
                       [1, 10],                 # tau_Q
                       [1, 10],                 # tau_U
                       [0.05, 0.5],             # d_SCSF
                       [0.007, 0.9],            # d_S
                       [0.7, 0.99],             # d_Q
                       [0.25, 0.85],            # d_U
                       [50, 5000],              # N_half
                       [1, 10],                 # S_PH
                       [1, 10],                 # S_PS
                       [5, 50],                 # S_PQ
                       [5, 50],                 # S_AU
                       [1, 10],                 # S_AH
                       [1, 10],                 # S_AS
                       [3, 30],                 # S_AM
                       [1000, 20000],           # S_SCSF
                       [1, 10],                 # S_KD
                       [1, 10],                 # k_sn
                       [2, 20],                 # k_nq
                       [1, 10],                 # k_nm
                       [1, 10],                 # k_ns
                       [1, 100],                # R_KU
                       [0.4, 0.95],             # I_crit
                       [2000, 40_000],          # K_crit
                       [0.1, 1],                # psi
                       [0.5, 1],                # d_M
                       [0.2, 0.9],              # d_MF
                       [3, 30],                 # S_KMD
                       [5, 50],                 # S_KQ
                       [1, 10],                 # C_QM
                       [1, 10],                 # C_MDM
                       [1, 10],                 # C_UM
                       [1000, 10_000],          # S_MF
                       [0.05, 0.95],            # omega
                       [1, 10],                 # C_UP
                       [0.1, 1],                # alpha
                       [0.1, 0.9],               # H_crit,
                       [30_000, 70_000]        # path_size
                       ]                
        }

    param_values_IHD = morris_sample.sample(problem_IHD, N = 500, num_levels = 8)       # number of samples generated = N * (P + 1) where P = # of parameters   
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

    output_to_analyze = np.zeros(param_values_IHD.shape[0])
    # ^ Doesn't store all outputs for memory efficiency purposes

    init_state = [

        30000,  # Quiescent HSPCs
        0,      # Proliferating HSPCs
        0,      # PAMPs (Pathogens)
        100,      # Pro-inflammatory Cytokines
        100,      # Anti-inflammatory Cytokines
        10000,  # Stem Cell Supporting Factors
        0,      # DAMPs (Tissue Damage)
        100,      # Activated leukocytes
        100,   # Stable leukocytes
        100,       # Suppressor leukocytes
        1,      # MDSC
        1000,  # MF

    ]       # init_state is mostly arbitary, need to check how this affects interpretability of SA results

    ext_stim = np.zeros((12, int(nTimesteps/0.01)))

    ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

    param_values_dict = [{**dict(zip(param_names[:-1], row[:-1])), 'N_oo': 2*10**7, 'k': 3, 'beta_N': 10**-3, 'd_P': 0.99, 'd_A': 0.99, 'g_N': 0.1, 'gamma': 0.0005, 'delta': 0.2} for row in param_values_IHD]           # convert generated samples array into an array of dictionaries

    print("! Computations starting now !")
    start_time = time.time()

    for i, X in enumerate(param_values_dict):

        ext_stim[2, 999] = param_values_IHD[i][-1]

        data = lin_sim(M3_beta.beta_model_3, X, init_state, nTimesteps, 0.01, ext_stim, ext_stim_m)   # simulate the digital system using sample input
        outputs = data[0]
        derivatives = data[1]

        chronic_arr = np.array(chronic_infection_flagger(1000, 200_000, 6000, 0.01, outputs[2], derivatives[2]))
        chronic_arr = chronic_arr[chronic_arr == 1]
        duration_of_chronic = len(chronic_arr) * 0.01

        output_to_analyze[i] = duration_of_chronic

        #IHD_out[i] = output[:, cols]
        print(f"Output {i+1} of {param_values_IHD.shape[0]} successfully calculated!")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Outputs generated for all input samples. Elapsed time: ", elapsed_time)

    print("Beginning calculation of morris indices for each timestep now...")    # sanity check
    start_time = time.time()


    SI_N = morris_analyze.analyze(problem_IHD, param_values_IHD, output_to_analyze, num_levels=8)
    df_SI_N = SI_N.to_df()

    df_SI_N['mu'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'chronic_duration_out', 'mu', f'mu_chronic_duration_{exp_num}.csv'), sep='\t')
    df_SI_N[['mu_star', 'mu_star_conf']].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'chronic_duration_out', 'mu_star', f'mu_star_chronic_duration_{exp_num}.csv'), sep='\t')
    df_SI_N['sigma'].to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'chronic_duration_out', 'sigma', f'sigma_chronic_duration_{exp_num}.csv'), sep='\t')


    end_time = time.time()
    print("All morris indices successfully calculate. Elapsed time: " + str(end_time - start_time))

# == Below code is used to plot mu_star against sigma =======
if generate_mu_vs_sigma:

    Filepath = os.path.join(script_dir, f'Experiment_{exp_num}')

    path = os.path.join(Filepath, 'mu_star_vs_sigma')
    if not os.path.exists(path):
        os.makedirs(path)


    for i, output in enumerate(output_names):
        
        df_mu_star = pd.read_csv(os.path.join(Filepath, f'{output}_out', 'mu_star', f'mu_star_{output}_{exp_num}.csv'), delimiter='\t')
        df_sigma = pd.read_csv(os.path.join(Filepath, f'{output}_out', 'sigma', f'sigma_{output}_{exp_num}.csv'), delimiter='\t')

        df_combined = pd.concat([df_mu_star, df_sigma], axis=1).iloc[:, [0, 1, 4]]\
        
        # this code here will sort and split the data into three to increase readability of graphs
        df_sorted = df_combined.sort_values(by='mu_star', ascending=True)
        df_sorted.to_csv(os.path.join(path, f'{output}_{exp_num}.csv'))

        tripoint = len(df_combined) // 3

        df_sorted_lowest = df_sorted.iloc[:tripoint]
        df_sorted_mid = df_sorted.iloc[tripoint:2*tripoint]
        df_sorted_highest = df_sorted.iloc[2*tripoint:]

        #print(df_combined.columns)             # sanity check
        

        for n, df in enumerate([df_sorted_lowest, df_sorted_mid, df_sorted_highest]):        # creates 2 separate figures

            plt.scatter(df["mu_star"], df["sigma"], color='blue')
            
            x_max = df["mu_star"].max()
            x_min = df['mu_star'].min()
            y_max = df["sigma"].max()
            y_min = df["sigma"].min()

            plt.xscale('log', base=10)
            plt.yscale('log', base=10)

            plt.plot(np.arange(0, x_max, 1), np.arange(0, x_max, 1), 'g--', label = '$f(x)=x$')         # plots f(x) = x
            plt.plot(np.arange(0, x_max, 1), np.arange(0, 2*x_max, 2), 'r--', label = '$f(x)=2x$')       # plots f(x) = 2x

            for j, row in df.iterrows():
                plt.annotate(row["Unnamed: 0"], (row["mu_star"], row["sigma"]), fontsize=8, xytext=(5,5), textcoords="offset points")

            plt.xlim(max(x_min, 0.01), x_max+10)
            plt.ylim(max(y_min, 0.01), y_max+100)

            plt.xlabel(r'$\mu^{*}$ (Mean Absolute Effect)')
            plt.ylabel(r'$\sigma$ (Variance)')
            plt.title(fr"$\mu^*$ vs $\sigma$ (${output_names_laTex[i]}$)")
            plt.legend(loc="upper left")
            plt.grid(True)

            plt.savefig(os.path.join(path, f'{output}_{exp_num}_{n}.png'))
            plt.close()
            #plt.show()
