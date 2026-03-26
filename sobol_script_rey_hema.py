import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from project_library import lin_sim_scipy
from project_library import csv_to_figure
from project_library import merge_figures_grid
from reynolds_project_func import reynolds_hema_ODE
import time
import pandas
import os

output_names = ['P', 'N', 'D', 'CA', 'M', 'I']

output_names_laTex = [r'P', r'N', r'D', r'C_{A}', r'M', r'I']

parameter_names = ['k_pm', 'k_mp', 's_m', 'mu_m', 'k_pg', 
                   'k_pn', 'k_np', 'k_nn', 's_i',
                   'mu_nr', 'mu_n', 'k_nd', 'k_dn', 'x_dn',
                   'mu_d', 'C_infty', 's_c', 'k_cn', 'k_cnd',
                   'mu_c', 'alpha_M', 'alpha_C', 'alpha_D', 'mu_md',
                   'alpha_MS1', 's_in']

parameter_names_laTex = [r'k_{pm}', r'k_{mp}', r's_{m}', r'\mu_{m}', r'k_{pg}',
                         r'k_{pn}', r'k_{np}', r'k_{nn}', r's_{i}',
                         r'\mu_{nr}', r'\mu_{n}', r'k_{nd}', r'k_{dn}', r'x_{dn}',
                         r'\mu_d', r'C_{\infty}', r's_{c}', r'k_{cn}', r'k_{cnd}',
                         r'\mu_{c}', r'\alpha_{M}', r'\alpha_{C}', r'\alpha_{D}', r'\mu_{md}',
                         r'\alpha_{MS1}', r's_{in}'] 

exp_num = 1             # trial number
generate_individual_figs = False
generate_merged_figs = True
run_sensitivity_analysis = False
nTimesteps = 50            # how many timesteps to run simulation
init_time = 25           # initial time to begin calculating sobol indices, NOT initial time of model
delta_t = 1             # timestep size for sobol index calculation, NOT the delta used in solver

# code directly below checks to see if the relevant folders already exist and if not creates them
script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hema_rey_SA_sobol')
if not os.path.exists(script_dir):
            os.makedirs(script_dir)

for out_name in output_names:
        
    for order in ['first', 'second', 'total']:

        filepath = os.path.join(script_dir, f'Experiment_{exp_num}', f'{out_name}_out', f'{order}_order')
        if not os.path.exists(filepath):
            os.makedirs(filepath)

if run_sensitivity_analysis:    # run sensitivity analysis and save relevant sensitivity index to .csv files with delimiter='\t\

    
    problem_IHD = {
            'num_vars': 26,
            'names': ['k_pm', 'k_mp', 's_m', 'mu_m', 'k_pg', 
                    'k_pn', 'k_np', 'k_nn', 's_nr',
                   'mu_nr', 'mu_n', 'k_nd', 'k_dn', 'x_dn',
                   'mu_d', 'C_infty', 's_c', 'k_cn', 'k_cnd',
                   'mu_c', 'alpha_M', 'alpha_C', 'alpha_D', 'mu_md',
                   'alpha_MS1', 's_in'],
            'bounds': [[1, 10],         # k_pm
                       [1, 10],         # k_mp
                       [1, 10],         # s_m
                       [1, 10],         # mu_m
                       [1, 10],         # k_pg
                       [1, 10],         # k_pn
                       [1, 10],         # k_np
                       [1, 10],         # k_nn
                       [1, 10],         # s_nr
                       [1, 10],         # mu_nr
                       [1, 10],         # mu_n
                       [1, 10],         # k_nd
                       [1, 10],         # k_dn
                       [1, 10],         # x_dn
                       [1, 10],         # mu_d
                       [1, 10],         # C_infty
                       [1, 10],         # S_c
                       [1, 10],         # k_cn
                       [1, 10],         # k_cnd
                       [1, 10],         # mu_c
                       [1, 10],         # alpha_M
                       [1, 10],         # alpha_C
                       [1, 10],         # alpha_D
                       [1, 10],         # mu_md
                       [1, 10],         # alpha_MS1
                       [1, 10]          # s_in
                       ]                
        }

    param_values_IHD = sobol_sample.sample(problem_IHD, 2048, calc_second_order=True)      # use n=2 when checking for errors, else n=2048 is a good standard
    print("Shape of the generated sample: ", param_values_IHD.shape)        # sanity check
    print("First few samples:")
    print(param_values_IHD[:3])

    IHD_out = np.zeros((param_values_IHD.shape[0], 6, int((nTimesteps-init_time))))     # this is where the necessary simulation outputs are stored for calculating indices
    # ^ Doesn't store all outputs for memory efficiency purposes
    IHD_out_test = []

    init_state = [

        0,
        0,
        0,
        0.125,
        0,
        0

    ]       # init_state is mostly arbitary, need to check how this affects interpretability of SA results

    param_values_dict = [{**dict(zip(parameter_names, row)), 'p_infty': 2*10**7} for row in param_values_IHD]           # convert generated samples array into an array of dicts, fix p_infty

    print("! Computations starting now !")
    nTick = 0
    start_time = time.time()

    cols = np.arange(init_time*100, nTimesteps*100, 100)
    print(cols)     # sanity check, make sure this is an array of the desired indices of the timesteps to save outputs from for index calculations
    t_eval = np.arange(0, nTimesteps, 0.01)
    

    for i, X in enumerate(param_values_dict):

        t, data = lin_sim_scipy(reynolds_hema_ODE, X, init_state, nTimesteps, 0.01, [0], [0], force_timesteps=True)   # simulate the digital system using sample input
        IHD_out[i] = data[:, cols]
        end_time = time.time()
        elapsed_time = end_time - start_time
        nTick = nTick + 1
        if nTick == 100:
             print(f"Output {i+1} of {param_values_IHD.shape[0]} successfully calculated. Elapsed time: {elapsed_time//60} min")
             nTick = 0

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Outputs generated for all input samples. Total elapsed time: ", elapsed_time)

    print("Beginning calculation of sobol indices for each timestep now...")    # sanity check
    start_time = time.time()

    for t in range(0, int(nTimesteps - init_time), 1):     # calculate sobol indices from init_time -> end of simulation with 1 timestep in between each

        SI_P = sobol_analyze.analyze(problem_IHD, IHD_out[:, 0, t], calc_second_order=True)    # int((init_time + t)/0.01)]
        SI_N = sobol_analyze.analyze(problem_IHD, IHD_out[:, 1, t], calc_second_order=True)
        SI_D = sobol_analyze.analyze(problem_IHD, IHD_out[:, 2, t], calc_second_order=True)
        SI_CA = sobol_analyze.analyze(problem_IHD, IHD_out[:, 3, t], calc_second_order=True)
        SI_M = sobol_analyze.analyze(problem_IHD, IHD_out[:, 4, t], calc_second_order=True)
        SI_I = sobol_analyze.analyze(problem_IHD, IHD_out[:, 5, t], calc_second_order=True)


        total_SI_P, first_SI_P, second_SI_P = SI_P.to_df()
        total_SI_N, first_SI_N, second_SI_N = SI_N.to_df()
        total_SI_D, first_SI_D, second_SI_D = SI_D.to_df()
        total_SI_CA, first_SI_CA, second_SI_CA = SI_CA.to_df()
        total_SI_M, first_SI_M, second_SI_M = SI_M.to_df()
        total_SI_I, first_SI_I, second_SI_I = SI_I.to_df()


        total_SI_P.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'total_order', f'total_SI_P_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        first_SI_P.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'first_order', f'first_SI_P_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        second_SI_P.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'P_out', 'second_order', f'second_SI_P_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        
        total_SI_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'total_order', f'total_SI_N_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        first_SI_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'first_order', f'first_SI_N_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        second_SI_N.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'N_out', 'second_order', f'second_SI_N_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        total_SI_D.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'D_out', 'total_order', f'total_Si_D_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        first_SI_D.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'D_out', 'first_order', f'first_Si_D_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        second_SI_D.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'D_out', 'second_order', f'second_Si_D_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        total_SI_CA.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'CA_out', 'total_order', f'total_Si_CA_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        first_SI_CA.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'CA_out', 'first_order', f'first_Si_CA_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        second_SI_CA.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'CA_out', 'second_order', f'second_Si_CA_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        total_SI_M.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'M_out', 'total_order', f'total_Si_M_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        first_SI_M.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'M_out', 'first_order', f'first_Si_M_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        second_SI_M.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'M_out', 'second_order', f'second_Si_M_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        total_SI_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'total_order', f'total_Si_I_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        first_SI_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'first_order', f'first_Si_I_{exp_num}_{int(t + init_time)}.csv'), sep='\t')
        second_SI_I.to_csv(os.path.join(script_dir, f'Experiment_{exp_num}', 'I_out', 'second_order', f'second_Si_I_{exp_num}_{int(t + init_time)}.csv'), sep='\t')

        print(f'Indices for timestep {int(t + init_time)} successfully calculated!')

    end_time = time.time()
    print("All sobol indices successfully calculate. Elapsed time: " + str(end_time - start_time))

if generate_individual_figs:    # generate time-series SI graphs for each permutation of (output, parameter)

    print("Figure generation starting now !")
    start_time = time.time()
    orders = ['first', 'total']

    for order in orders:

        filepath = os.path.join(script_dir, f'Experiment_{exp_num}')      # this path should be the branch containing all the output folders which themselves contain the relevant .csv's
        csv_to_figure(output_names, output_names_laTex, parameter_names, parameter_names_laTex, nTimesteps - init_time, init_time, 2, order, exp_num, filepath)      # this function takes care of placing the figures in the correct output folders automatically

        print(f"Figures succesfully generated for {order}-order indices.")

    
    end_time = time.time()
    print("Individual time-series figures generated for each output and parameter. Elapsed time: " + str(end_time-start_time))

if generate_merged_figs:        # merge time-series SI graphs for each output by parameter into one .pngS

    print("Merged figures generation starting now !")
    start_time = time.time()
    for order in ['first', 'total']:

        merge_figures_grid(2, 3, 800, 600, exp_num, order, output_names, parameter_names, directory=script_dir)
        print(f"Merged figures successfully generated for {order}-order indices.")
    
    end_time = time.time()
    print("Merged figures successfully generated. Elapsed time: " + str(end_time-start_time))
