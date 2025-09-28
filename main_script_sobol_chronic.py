import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from project_library import lin_sim
from project_library import csv_to_figure
from project_library import merge_figures_grid
from project_library import chronic_infection_flagger
import time
import pandas
from pathlib import Path
import M3_beta

# output_names will be used for formatting file names / output_names_laTex will be used when generating matplotlib plots
# same rule applies to param_names / param_names_laTex
# see if-else statements below to determine which preset to use

output_names = ['chronic_duration']

output_names_laTex = ['Chronic Duration']

param_names = ['H_crit', 'K_crit', 'S_AH', 'S_PH', 'S_SCSF', 'alpha', 'dH', 'd_SCSF', 'd_S', 'k_H', 'k_ns', 'psi', 'theta_N', 'path_size']
    
param_names_laTex = ['H_{crit}', 'K_{crit}', 'S_{AH}', 'S_{PH}', 'S_{SCSF}', '\\alpha', 'd_{H}', 'd_{SCSF}', 'd_{S}', 'k_{H}', 'k_{ns}', '\psi', '\\theta_{N}', 'path_size']



exp_num = 'sobol_chronic_final_debug'             # trial number, doesn't have to be number

nTimesteps = 400            # how many timesteps to run simulation
t_infection = 150           # timestep to introduce infection

problem_IHD = {
        'num_vars': 14,
        'names': param_names,
        'bounds': [
                [0.1, 0.9],    # H_crit
                [20_000, 200_000],    # K_crit
                [1, 10],           # S_AH
                [1, 10],           # S_PH
                [1000, 20_000],    # S_SCSF
                [0.1, 1],          # alpha
                [0.05, 0.8],       # dH
                [0.05, 0.5],       # dSCSF
                [0.007, 0.9],      # d_S
                [1, 10],           # k_H
                [1, 10],           # k_ns
                [0.1, 1],          # psi
                [5000, 50_000],     # theta_N
                [30_000, 70_000]   # path_size
        ]
    }

param_values_IHD = sobol_sample.sample(problem_IHD, N=2048, calc_second_order=True)      # use n=2 when checking for errors, else use n=2048; number of samples = n * (2(num of params) + 2)
print("Shape of the generated sample list: ", param_values_IHD.shape)
'''print("First few samples:")
print(param_values_IHD[:3])         # sanity check'''

IHD_out = np.zeros(param_values_IHD.shape[0])       # will contain the chronic infection durations associated with each run

init_state = [

    25000,  # Quiescent HSPCs
    5000,      # Proliferating HSPCs
    1,      # PAMPs (Pathogens)
    100,      # Pro-inflammatory Cytokines
    100,      # Anti-inflammatory Cytokines
    10000,  # Stem Cell Supporting Factors
    1,      # DAMPs (Tissue Damage)
    100,      # Activated leukocytes
    100,   # Stable leukocytes
    100,       # Suppressor leukocytes
    1,      # MDSC
    1000,  # MF

]

ext_stim = np.zeros((12, int(nTimesteps/0.01)))
ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

fixed_parameters = {

        'theta_K' : 10_000,
        'tau_Q' : 1,
        'tau_U' : 1,    
        'd_Q' : 0.95,
        'd_U' : 0.5,
        'd_P' : 0.99,
        'd_A' : 0.99,
        'g_N' : 0.10,
        'N_oo' : 2 * 10**7,
        'N_half' : 2500,
        'S_PS' : 1,
        'S_PQ' : 10,
        'S_AU' : 15,
        'S_AS' : 0,
        'S_AM' : 8,
        'S_KD' : 1,
        'k_sn' : 3,
        'k_nq' : 20,
        'k_nm' : 3,
        'R_KU' : 10,
        'I_crit' : 0.4,
        'k' : 3,
        'd_M' : 9/10,
        'd_MF' : 0.3,
        'S_KMD' : 2,
        'S_KQ' : 4,
        'C_QM' : 1,
        'C_MDM' : 1,
        'C_UM' : 1/3,
        'S_MF' : 500,
        'omega' : 0.5,
        'C_UP' : 2,
        'beta_N': 10**-3,
        'gamma': 0.0005,
        'delta': 0.2,

    }

param_values_dict = [{**dict(zip(param_names[:-1], row[:-1])), **fixed_parameters} for row in param_values_IHD]           # convert generated samples array into an array of dicts
pathogen_input = [row[-1] for row in param_values_IHD]

print("! Computations starting now !")
start_time = time.time()


for i, X in enumerate(param_values_dict):

    start_time_for_one = time.time()

    stim = ext_stim
    stim[2, int(t_infection/0.01)] = pathogen_input[i]

    '''if i < 5:   # sanity check
            print(stim[2, int(t_infection/0.01)-1 : int(t_infection/0.01)+2])'''

    data = lin_sim(M3_beta.beta_model_3, X, init_state, nTimesteps, 0.01, stim, ext_stim_m)   # simulate the digital system using sample input; index 0 -> model output, index 1 -> derivative output

    outputs = data[0]
    derivatives = data[1]

    chronic_arr = np.array(chronic_infection_flagger(1000, 200_000, 6000, 0.01, outputs[2], derivatives[2]))
    chronic_arr = chronic_arr[chronic_arr == 1]
    duration_of_chronic = len(chronic_arr) * 0.01

    IHD_out[i] = duration_of_chronic

    end_time_for_one = time.time()
    print(f"Output {i+1} of {param_values_IHD.shape[0]} successfully calculated!")

    if i == 0:
            print(f"Time for one output: {end_time_for_one - start_time_for_one} sec")
            print(f"Approximate total time to completion: {(param_values_IHD.shape[0] * (end_time_for_one - start_time_for_one)) / 3600} hrs")

end_time = time.time()
elapsed_time = end_time - start_time
print("Outputs generated for all input samples. Elapsed time: ", elapsed_time)


# code directly below checks to see if the relevant folders already exist and if not creates them
path = Path.cwd() / 'Sobol_Analysis_Final'
if not Path.exists(path):
    Path.mkdir(path)


# code below is WIP code for sobol on chronic duration

total_SI, first_SI, second_SI = sobol_analyze.analyze(problem_IHD, IHD_out, calc_second_order=True).to_df()

total_SI.to_csv(path / 'chronic_total_SI.csv', index=True)
first_SI.to_csv(path / 'chronic_first_SI.csv', index=True)
second_SI.to_csv(path / 'chronic_second_SI.csv', index=True)

print("All sobol indices successfully calculate")


