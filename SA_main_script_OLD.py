import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from SALib.test_functions import Ishigami
from project_library import linear_sim
from project_library import csv_to_figure
from project_library import merge_figures_grid
import time
import pandas
import os

output_names = ['H', 'S', 'I', 'N', 'Q', 'U']

param_names = ['gn', 'knq', 'kns', 'ktn', 'w', 'pcrit', 'spq', 'sph',
                'sps', 'sas', 'sah', 'sau', 'Is', 'As', 'Ic', 'Ac',
                'y', 'ds', 'dp', 'da', 'dq', 'du', 'H_init', 'N_init']  # important: parameter names MUST be in same order as they appear in the .csv file

exp_num = 3
generate_individual_figs = True
generate_merged_figs = True
run_sensitivity_analysis = True
nTimesteps = 100
offset = 0

# code block directly below is using linear_sim function to evaluate sensitivity on the IHD model 
if run_sensitivity_analysis:

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

    # init_values = [1000, 0, 4000, 0, 0.9, 0.1, 1800, 300, 1] # chosen mostly arbitrarily

    param_values_IHD = sobol_sample.sample(problem_IHD, 2)      # used n=2 when checking for errors, else use either n=1024 or n=2048
    print("Shape of the generated sample: ", param_values_IHD.shape)
    print("First few samples:")
    print(param_values_IHD[:3])

    init_values = np.hstack((param_values_IHD[:, -2:-1], np.tile(np.array([0, 4000, 0, 0.9, 0.1, 1800, 300, 1]), (param_values_IHD.shape[0], 1))))
    print("First few initial value sets: ")         # sanity check
    print(init_values[3,:])
    print("Shape of init_value array: ", init_values.shape)

    IHD_out = np.zeros((param_values_IHD.shape[0], 9, nTimesteps))      # 100 because want SA results for entire 100 hours of simulation

    print("! Computations starting now !")
    start_time = time.time()
    for i, X in enumerate(param_values_IHD):

        output = linear_sim(init_values[i], X[:23], 0.01, 100, [50], [X[-1]])   # simulate the digital system using sample input

        for x in range(nTimesteps-offset):

            IHD_out[i, :, x] = output[:, int(x/0.01)]

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Outputs generated for all input samples. Elapsed time: ", elapsed_time)

    print("Beginning calculation of sobol indices for each timestep now...")    # sanity check
    start_time = time.time()

    # code directly below checks to see if the relevant folders already exist and if not creates them
    script_dir = os.path.dirname(os.path.abspath(__file__))

    for out_name in output_names:
            
        for order in ['first', 'second', 'total']:

            filepath = os.path.join(script_dir, f'Experiment_{exp_num}', f'{out_name}_out', f'{order}_order')
            if not os.path.exists(filepath):
                os.makedirs(filepath)


    for t in range(nTimesteps):     # must update to account for possible offset

        SI_H = sobol_analyze.analyze(problem_IHD, IHD_out[:, 0, t])
        SI_N = sobol_analyze.analyze(problem_IHD, IHD_out[:, 1, t])
        SI_S = sobol_analyze.analyze(problem_IHD, IHD_out[:, 4, t])
        SI_Q = sobol_analyze.analyze(problem_IHD, IHD_out[:, 5, t])
        SI_U = sobol_analyze.analyze(problem_IHD, IHD_out[:, 6, t])
        SI_I = sobol_analyze.analyze(problem_IHD, IHD_out[:, 8, t])

        total_Si_H, first_Si_H, second_Si_H = SI_H.to_df()
        total_Si_N, first_Si_N, second_Si_N = SI_N.to_df()
        total_Si_S, first_Si_S, second_Si_S = SI_S.to_df()
        total_Si_Q, first_Si_Q, second_Si_Q = SI_Q.to_df()
        total_Si_U, first_Si_U, second_Si_U = SI_U.to_df()
        total_Si_I, first_Si_I, second_Si_I = SI_I.to_df()

        total_Si_H.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'H_out', 'total_order', f'total_Si_H_{exp_num}_{t}.csv'), sep='\t')
        first_Si_H.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'H_out', 'first_order', f'first_Si_H_{exp_num}_{t}.csv'), sep='\t')
        second_Si_H.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'H_out', 'second_order', f'second_Si_H_{exp_num}_{t}.csv'), sep='\t')

        total_Si_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'total_order', f'total_Si_N_{exp_num}_{t}.csv'), sep='\t')
        first_Si_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'first_order', f'first_Si_N_{exp_num}_{t}.csv'), sep='\t')
        second_Si_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'second_order', f'second_Si_N_{exp_num}_{t}.csv'), sep='\t')

        total_Si_S.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'total_order', f'total_Si_S_{exp_num}_{t}.csv'), sep='\t')
        first_Si_S.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'first_order', f'first_Si_S_{exp_num}_{t}.csv'), sep='\t')
        second_Si_S.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'S_out', 'second_order', f'second_Si_S_{exp_num}_{t}.csv'), sep='\t')

        total_Si_Q.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'total_order', f'total_Si_Q_{exp_num}_{t}.csv'), sep='\t')
        first_Si_Q.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'first_order', f'first_Si_Q_{exp_num}_{t}.csv'), sep='\t')
        second_Si_Q.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'Q_out', 'second_order', f'second_Si_Q_{exp_num}_{t}.csv'), sep='\t')

        total_Si_U.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'total_order', f'total_Si_U_{exp_num}_{t}.csv'), sep='\t')
        first_Si_U.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'first_order', f'first_Si_U_{exp_num}_{t}.csv'), sep='\t')
        second_Si_U.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'U_out', 'second_order', f'second_Si_U_{exp_num}_{t}.csv'), sep='\t')

        total_Si_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'total_order', f'total_Si_I_{exp_num}_{t}.csv'), sep='\t')
        first_Si_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'first_order', f'first_Si_I_{exp_num}_{t}.csv'), sep='\t')
        second_Si_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'second_order', f'second_Si_I_{exp_num}_{t}.csv'), sep='\t')

    end_time = time.time()
    print("All sobol indices successfully calculate. Elapsed time: " + str(end_time - start_time))

if generate_individual_figs:

    start_time = time.time()
    orders = ['first', 'total']

    for order in orders:

        filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'Experiment_{exp_num}')      # this path should be the branch containing all the folders which contain the relevant .csv's

        csv_to_figure(output_names, param_names, nTimesteps, offset, 2, order, exp_num, filepath)      # this function takes care of placing the figures in the correct output folders automatically
        # init_time = offset
    
    end_time = time.time()
    print("Individual time-series figures generated for each output and parameter. Elapsed time: " + str(end_time-start_time))

if generate_merged_figs:

    start_time = time.time()
    for order in ['first', 'total']:

        merge_figures_grid(2, 3, 800, 600, 3, order, output_names, param_names)
    
    end_time = time.time()
    print("Merged figures successfully generated. Elapsed time: " + str(end_time-start_time))\
    
# note that second_order figures must be handled separately due to the difference in data shape
