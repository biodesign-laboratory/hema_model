"""
csv is the default I/O format.

This file is like a debugging tool and a means to quickly explore some outputs
before writing the main plotting code in generate_figures.py.

This file allows the use of terminal flags, which are not only more
user-friendly, but directly compatible with running simulations on
the cluster.

I'll work on getting this file to act like the original sim_V2.py.
"""

import matplotlib.pyplot as plt
import numpy as np
import project_library as PL
import time
from M3_beta import beta_model_3
import pandas as pd
from pathlib import Path
import csv

import argparse
import os



# checks to see if the preset folder exists and if not makes it
PATH_PRESETS = Path.cwd() / 'presets'
if not PATH_PRESETS.exists():
    PATH_PRESETS.mkdir(parents=True, exist_ok=True)


def load_hyper_parameters_(read_from_hyper_loc):
    with open(read_from_hyper_loc, mode="r", encoding="utf-8") as file:

        reader = csv.reader(file)

        #hyper_from_csv = {rows[0]: float(rows[1]) for rows in reader}

        hyper_from_csv = {}
        for row in reader:
            key = row[0]
            value = row[1]

            # the keys delayed_infection, HSPC_boost, and SCSF_boost are bools
            if value.lower() == "true":
                value = True
            elif value.lower() == "false":
                value = False
            else:
                if key == 'runs':
                    # 'runs' must be int
                    value = int(value)
                else:
                    try:
                        value = float(value)
                    except ValueError:
                        pass

            hyper_from_csv[key] = value

    print("Hyperparameters successfully read from .csv provided")
    # print(hyper_from_csv)   # sanity check

    return hyper_from_csv

def save_hyper_parameters_(hyper_fname,hyper_parameters):

    with open(hyper_fname, "w", newline='') as file:
        writer = csv.writer(file)
        
        for key, value in hyper_parameters.items():
            writer.writerow([key, value])

    print("Hyperparameters successfully saved to .csv")


def load_init_states_(read_from_init_loc):
    
    with open(read_from_init_loc, mode="r", encoding="utf-8") as file:

        reader = csv.reader(file)
        init_from_csv = {rows[0]: float(rows[1]) for rows in reader}

    print("Initial values successfully read from .csv")

    return init_from_csv


def save_init_states_(init_fname,init_dict):

    with open(init_fname, "w", newline='') as file:
        writer = csv.writer(file)
        
        for key, value in init_dict.items():
            writer.writerow([key, value])
            
        print("Initial values successfully saved to .csv")


def load_params_(read_from_param_loc):
    with open(read_from_param_loc, mode="r", encoding="utf-8") as file:

        reader = csv.reader(file)
        parameters = {rows[0]: float(rows[1]) for rows in reader}

    
    print("Model parameters successfully read from provided .csv")

    return parameters

def save_params_(param_fname,parameters):
    
    with open(param_fname, "w", newline='') as file:
        writer = csv.writer(file)

        for key, value in parameters.items():
            writer.writerow([key, value])

    print("Model parameters successfully saved to .csv")


def get_default_hyp():
    default_hyp = {
        'delta_t' : 0.01,
        't_final' : 400,
        'delayed_infection' : False,
        't_infection_start' : 100,
        't_infection_end' : 200,
        'path_size_default' : 40_000,
        'nosocomial_size' : 10_000,
        'nosocomial_start' : 250,
        'nosocomial_end' : 350,
        'path_increment' : 2000,
        'HSPC_boost' : False,
        'HSPC_default' : 0,
        'HSPC_increment' : 5_000,
        'HSPC_boost_time' : 110,
        'SCSF_boost' : False,
        'SCSF_default' : 0,
        'SCSF_increment' : 0,
        'SCSF_boost_time' : 110
    }

    return default_hyp

def get_default_inits():
    default_inits = {
        'HQ':10_000,
        'HM':5_000,
        'N':0,
        'P':10_000,
        'A':10_000,
        'SCSF':10_000,
        'K':1,
        'Q':500,
        'S':500,
        'U':500,
        'MDSC':1,
        'MF':5_000
    }
    return default_inits


def get_default_params():
    default_params = {
        'k_H' : 3,
        'dH' : 0.05,
        'theta_N' : 100_000,
        'theta_K' : 50_000,
        'tau_Q' : 1,
        'tau_U' : 2,
        'd_SCSF' : 0.3,
        'd_S' : 0.05,
        'd_Q' : 0.95,
        'd_U' : 0.5,
        'd_P' : 0.95,
        'd_A' : 0.95,
        'g_N' : 0.10,
        'N_oo' : 2 * 10**7,
        'N_half' : 2500,
        'S_PH' : 2,
        'S_PS' : 5,
        'S_PQ' : 10,
        'S_AU' : 15,
        'S_AH' : 0,
        'S_AS' : 0,
        'S_AM' : 12,
        'S_SCSF' : 10000,
        'S_KD' : 1,
        'k_sn' : 3,
        'k_nq' : 10,
        'k_nm' : 3,
        'k_ns' : 0.5,
        'R_KU' : 10,
        'I_crit' : 0.8,
        'K_crit' : 150_000,
        'k' : 1,
        'psi' : 1,
        'd_M' : 9/10,
        'd_MF' : 0.3,
        'S_KMD' : 2,
        'S_KQ' : 4,
        'C_QM' : 1,
        'C_MDM' : 1,
        'C_UM' : 1/3,
        'S_MF' : 1000,
        'omega' : 0.8,
        'C_UP' : 2,
        'alpha': 1/3,
        'beta_N': 10**-3
    }
    return default_params


def get_parameters(read_from_param_loc,read_from_param,run_number,save_parameters):
    if not(os.path.isfile(read_from_param_loc)) or not(read_from_param):
        # set model parameters manually

        parameters = get_default_params()

        if save_parameters:
            param_fname = PATH_PRESETS / f'parameter_preset_{run_number}.csv'
            save_params_(param_fname,parameters)

    else:           # set model parameters by reading from .csv
        # Warning: if .csv key-values do not match selected model format, error in sim will occur
        parameters = load_params_(read_from_param_loc)

    return parameters


def get_init_states(read_from_init_loc,read_from_init,run_number,save_init_states):
    if not(os.path.isfile(read_from_init_loc)) or not(read_from_init):
        # set initial values manually
        # this dict will be saved as a .csv if save_init_states = True
        init_dict = get_default_inits()

        if save_init_states:        # save initial values to a .csv preset
            init_fname = PATH_PRESETS / f'init_val_preset_{run_number}.csv'
            save_init_states_(init_fname,init_dict)

    else:               # set initial values by reading from .csv

        # Warning: if .csv key-values do not match selected model format, error in sim will occur
        init_dict = load_init_states_(read_from_init_loc)

    return init_dict


def get_hyp(read_from_hyper_loc,read_from_hyper,run_number,save_hyperparams):
    if not(os.path.isfile(read_from_hyper_loc)) or not(read_from_hyper):
        hyp = get_default_hyp()
        if save_hyperparams:
            # file name and location to save to
            hyper_fname = PATH_PRESETS / f'hyper_preset_{run_number}.csv'
            save_hyper_parameters_(hyper_fname,hyp)

    else:           # set hyperparameters by reading from a .csv
        hyp = load_hyper_parameters_(read_from_hyper_loc)

    return hyp


def main():
    
    parser = argparse.ArgumentParser(description='Run simulations')

    parser.add_argument('--save-init-states', action=argparse.BooleanOptionalAction,
                        help='Save initial conditions.',default=True)

    parser.add_argument('--save-parameters', action=argparse.BooleanOptionalAction,
                        help='Save parameters.',default=True)
    
    parser.add_argument('--save-hyper', action=argparse.BooleanOptionalAction,
                        help='Save hyper parameters.',default=True)

    parser.add_argument('--read-hyper', action='store_true',help='Read hyper parameters.')
    parser.add_argument('--read-param', action='store_true',help='Read parameters.')
    parser.add_argument('--read-init', action='store_true',help='Read initial conditions.')

    parser.add_argument('--save-solution', action='store_true',help='Save outputs to csv file.')

    args = parser.parse_args()
    
    run_number = 'chronic_1'                # used in file names, doesn't have to be a number
    output_to_csv = args.save_solution

    save_init_states = args.save_init_states # save initial conditions to .csv
    save_parameters = args.save_parameters   # save parameters to .csv
    save_hyperparams = args.save_hyper       # save hyperparameters to .csv

    read_from_hyper = args.read_hyper # whether hyperparameters should be loaded
    read_from_param = args.read_param # whether model parameters should be loaded
    read_from_init = args.read_init  # whether initial values should be loaded

    print(args)
    print(save_init_states)

    #default_hyp = get_default_hyp()
    #default_inits = get_default_inits()    
    #default_params = get_default_params() 
    

    # .csv containing hyperparameters to read from, put the path here
    read_from_hyper_loc = Path.cwd() / 'presets' / 'hyper_preset_{}.csv'.format(run_number)
    # .csv containing model paramaters to read from, put the path here
    read_from_param_loc = Path.cwd() / 'presets' / 'parameter_preset_{}.csv'.format(run_number)
    # .csv containing model initial values to read from, put the path here
    read_from_init_loc = Path.cwd() / 'presets' / 'init_val_preset_{}.csv'.format(run_number)

    # the following arguments are used to input artifical quiescent HSPCs, not included in csv
    #################################### hyperparameters here    


    hyp = get_hyp(read_from_hyper_loc,read_from_hyper,run_number,save_hyperparams)
    
    delta_t = hyp['delta_t']
    t_final = hyp['t_final']
    delayed_infection = hyp['delayed_infection']
    t_infection_start = hyp['t_infection_start']
    t_infection_end = hyp['t_infection_end']
    path_size_default = hyp['path_size_default']
    nosocomial_size = hyp['nosocomial_size']
    nosocomial_start = hyp['nosocomial_start']
    nosocomial_end = hyp['nosocomial_end']
    path_increment = hyp['path_increment']

    HSPC_boost = hyp['HSPC_boost']
    HSPC_default = hyp['HSPC_default']
    HSPC_increment = hyp['HSPC_increment']
    HSPC_boost_time = hyp['HSPC_boost_time']

    SCSF_boost = hyp['SCSF_boost']
    SCSF_default = hyp['SCSF_default']
    SCSF_increment = hyp['SCSF_increment']
    SCSF_boost_time = hyp['SCSF_boost_time']

    #timesteps = np.arange(stop=t_final, step=delta_t)

    #################################### initial states here

    init_dict = get_init_states(read_from_init_loc,read_from_init,run_number,save_init_states)

    init_state = [
        init_dict['HQ'],
        init_dict['HM'],
        init_dict['N'],
        init_dict['P'],
        init_dict['A'],
        init_dict['SCSF'],
        init_dict['K'],
        init_dict['Q'],
        init_dict['S'],
        init_dict['U'],
        init_dict['MDSC'],
        init_dict['MF']
    ]

    #################################### parameters here
    parameters = get_parameters(read_from_param_loc,read_from_param,run_number,save_parameters)


        
    # ----- 0. Generating arrays before running solver -------------
    
    num_outputs = 13
    #ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

    #ext_stimuli = np.zeros((num_outputs-1, len(timesteps)))

    #if delayed_infection == False:                          # one-time pathogen input
    #    ext_stimuli[2, int(t_infection_start/delta_t)] = path_size_default
    #    ext_stimuli[2, int(nosocomial_start/delta_t)] = nosocomial_size # nosocomial

    #else:                                                   # delayed input
    #    for j in np.arange(int(t_infection_start/delta_t), int(t_infection_end/delta_t)):

    #        ext_stimuli[2, j] = path_size_default / len(np.arange(int(t_infection_start/delta_t), int(t_infection_end/delta_t)))

    #    for j in np.arange(int(nosocomial_start/delta_t), int(nosocomial_end/delta_t)):

    #        ext_stimuli[2, j] = (nosocomial_size) / len(np.arange(int(nosocomial_start/delta_t), int(nosocomial_end/delta_t)))

    #if HSPC_boost:
    #    
    #    ext_stimuli[0, int(HSPC_boost_time/delta_t)] = HSPC_default + HSPC_increment*i

    #if SCSF_boost:
    #    ext_stimuli[5, int(SCSF_boost_time/delta_t)] = SCSF_default + SCSF_increment*i

    #outputs = np.zeros((num_outputs, len(timesteps)))


    # ------- 1. Run ODE solver -------

    stim_times = [100]
    stim_sizes = [2000]

    start = time.time()
    t,data = PL.lin_sim_scipy(beta_model_3, parameters, init_state, t_final,
                              delta_t, stim_times, stim_sizes)
    I = PL.calculate_I(data[3], data[4], data[6], data[2], parameters['theta_N'],
                       parameters['theta_K'], parameters['k'])

    outputs = np.concatenate([data,I.reshape(1,-1)])

    print(f"Run output successfully computed")
    end = time.time()
    print(f'Execution successful. Time elapsed: {end-start}s')

    # ---- 2. Save output to .csv's and plot the simulations -------
    path_r = Path.cwd() / 'Runs'

    if not Path.exists(path_r):
        Path.mkdir(path_r)
        
    path_m = path_r / 'Model 3'

    if not Path.exists(path_m):
        Path.mkdir(path_m)

    path_e = path_m / f'Exp_{run_number}'

    if not Path.exists(path_e):
        Path.mkdir(path_e)

    if output_to_csv:
        print("Saving outputs to .csv files now.")
        start = time.time()

        df = pd.DataFrame(outputs.T, columns=['HQ', 'HM', 'N', 'P', 'A', 'SCSF', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I'], index=list(map(str, t)))
        df.to_csv(path_e / f'SIM_{run_number}_output.csv')
        
        end = time.time()
        print(f"Outputs successfully saved to .csv's. Time elapsed: {end-start}s")

    fig1, axs = plt.subplots(3, 4, figsize=(8,7))
    fig2, axs2 = plt.subplots(1, 1, figsize=(8,2))

    axs1 = np.asarray(axs).flatten()

    # graph titles here
    titles = [
        '$H_Q(t)$',
        '$H_M(t)$',
        '$N(t)$',
        '$P(t)$',
        '$A(t)$',
        '$SCSF(t)$',
        '$K(t)$',
        '$Q(t)$',
        '$S(t)$',
        '$U(t)$',
        '$MDSC(t)$',
        '$MF(t)$', 
        '$I(t)$'
    ]


    for i,out in enumerate(outputs[:-1]):
        axs1[i].plot(t, out)
        axs1[i].set_title(titles[i])

    axs2.plot(I)
    
    plt.tight_layout()

    # ----------------------------------------------------------

    #fig1.savefig(path_e / f'sim_{run_number}_Hq_Hp_N.png', dpi=300)
    #fig2.savefig(path_e / f'sim_{run_number}_P_A_SCSF.png', dpi=300)
    #fig3.savefig(path_e / f'sim_{run_number}_K_Q_S.png', dpi=300)
    #fig4.savefig(path_e / f'sim_{run_number}_U_MDSC_MF.png', dpi=300)
    #fig5.savefig(path_e / f'sim_{run_number}_I.png', dpi=300)

    plt.show()


if __name__ == "__main__":
    main()
