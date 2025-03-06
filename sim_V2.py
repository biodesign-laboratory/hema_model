import matplotlib.pyplot as plt
import numpy as np
import project_library as PL
import time
from M3_beta import beta_model_3
import pandas as pd
from pathlib import Path
import csv

run_number = 'chronic_1'                # used in file names, doesn't have to be a number
model = 3                               # 2 - model 2 (previous, no MDSCs), 3 - model 3 (MDSCs)
output_to_csv = False     

save_init_states = True                # save initial conditions to .csv
save_parameters = True                # save parameters to .csv
save_hyperparams = True                # save hyperparameters to .csv

read_from_hyper = False              # whether hyperparameters should be loaded from .csv
read_from_param = False              # whether model parameters should be loaded from .csv
read_from_init = False               # whether initial values should be loaded from .csv

read_from_hyper_loc = Path.cwd() / 'presets' / 'hyper_preset_debug.csv'         # .csv containing hyperparameters to read from, put the path here
read_from_param_loc = Path.cwd() / 'presets' / 'parameter_preset_debug.csv'        # .csv containing model paramaters to read from, put the path here
read_from_init_loc = Path.cwd() / 'presets' / 'init_val_preset_debug.csv'         # .csv containing model initial values to read from, put the path here

# the following arguments are used to input artifical quiescent HSPCs, not included in csv



# checks to see if the preset folder exists and if not makes it
path = Path.cwd() / 'presets'
if not path.exists():
    path.mkdir(parents=True, exist_ok=True)

# hyperparameters here
if not read_from_hyper:             # set hyperparameters manually
    runs = 10
    delta_t = 0.01
    t_final = 400                       # 672 hours = 4 weeks (timescale is arbitrary but I still use it to gauge model behavior)
    delayed_infection = False            # use this to make pathogen input spread out over time period [t_infection_start, t_infection_end)
    t_infection_start = 100             
    t_infection_end = 200               # only used if delayed_infection = True
    path_size_default = 40_000
    nosocomial_size = 10000
    nosocomial_start = 250
    nosocomial_end = 350                # only used if delayed_infection = True
    path_increment = 2000

    HSPC_boost = False
    HSPC_default = 0
    HSPC_increment = 5000
    HSPC_boost_time = 110

    SCSF_boost = False
    SCSF_default = 0
    SCSF_increment = 0
    SCSF_boost_time = 110

    if save_hyperparams:
        
        hyper_parameters = {
            'runs' : runs,
            'delta_t' : delta_t,
            't_final' : t_final,
            'delayed_infection' : delayed_infection,
            't_infection_start' : t_infection_start,
            't_infection_end' : t_infection_end,
            'path_size_default' : path_size_default,
            'nosocomial_size' : nosocomial_size,
            'nosocomial_start' : nosocomial_start,
            'nosocomial_end' : nosocomial_end,
            'path_increment' : path_increment,
            'HSPC_boost' : HSPC_boost,
            'HSPC_default' : HSPC_default,
            'HSPC_increment' : HSPC_increment,
            'HSPC_boost_time' : HSPC_boost_time,
            'SCSF_boost' : SCSF_boost,
            'SCSF_default' : SCSF_default,
            'SCSF_increment' : SCSF_increment,
            'SCSF_boost_time' : SCSF_boost_time

        }

        hyper_fname = path / f'hyper_preset_{run_number}.csv'          # file name and location to save to

        with open(hyper_fname, "w", newline='') as file:
            writer = csv.writer(file)
            
            for key, value in hyper_parameters.items():
                writer.writerow([key, value])
        
        print("Hyperparameters successfully saved to .csv")

else:           # set hyperparameters by reading from a .csv

    with open(read_from_hyper_loc, mode="r", encoding="utf-8") as file:

        reader = csv.reader(file)

        #hyper_from_csv = {rows[0]: float(rows[1]) for rows in reader}

        hyper_from_csv = {}
        for row in reader:
            key = row[0]
            value = row[1]

            if value.lower() == "true":     # the keys delayed_infection, HSPC_boost, and SCSF_boost are bools
                value = True
            elif value.lower() == "false":
                value = False
            else:
                if key == 'runs':
                    value = int(value)      # 'runs' specifically needs to be an int
                else:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
            
            hyper_from_csv[key] = value

    runs = hyper_from_csv['runs']
    delta_t = hyper_from_csv['delta_t']
    t_final = hyper_from_csv['t_final']
    delayed_infection = hyper_from_csv['delayed_infection']
    t_infection_start = hyper_from_csv['t_infection_start']
    t_infection_end = hyper_from_csv['t_infection_end']
    path_size_default = hyper_from_csv['path_size_default']
    nosocomial_size = hyper_from_csv['nosocomial_size']
    nosocomial_start = hyper_from_csv['nosocomial_start']
    nosocomial_end = hyper_from_csv['nosocomial_end']
    path_increment = hyper_from_csv['path_increment']

    HSPC_boost = hyper_from_csv['HSPC_boost']
    HSPC_default = hyper_from_csv['HSPC_default']
    HSPC_increment = hyper_from_csv['HSPC_increment']
    HSPC_boost_time = hyper_from_csv['HSPC_boost_time']

    SCSF_boost = hyper_from_csv['SCSF_boost']
    SCSF_default = hyper_from_csv['SCSF_default']
    SCSF_increment = hyper_from_csv['SCSF_increment']
    SCSF_boost_time = hyper_from_csv['SCSF_boost_time']


    print("Hyperparameters successfully read from .csv provided")
    # print(hyper_from_csv)   # sanity check

timesteps = np.arange(stop=t_final, step=delta_t)

# initial states here
if not read_from_init:        # set initial values manually
    
    if model == 2:
        init_state = [

            10000,  # Quiescent HSPCs
            5000,      # Proliferating HSPCs
            0,      # PAMPs (Pathogens)
            10000,      # Pro-inflammatory Cytokines
            10000,      # Anti-inflammatory Cytokines
            10000,  # Stem Cell Supporting Factors
            0,      # DAMPs (Tissue Damage)
            0,      # Activated leukocytes
            600,   # Stable leukocytes
            0       # Suppressor leukocytes

        ]

        output_keys = ['HQ', 'HM', 'N', 'P', 'A', 'SCSF', 'K', 'Q', 'S', 'U']

        init_dict = dict(zip(output_keys, init_state))              # this dict will be saved as a .csv if save_init_states = True

    elif model == 3:
        init_state = [

            10000,  # Quiescent HSPCs
            5000,      # Proliferating HSPCs
            0,      # PAMPs (Pathogens)
            10000,      # Pro-inflammatory Cytokines
            10000,      # Anti-inflammatory Cytokines
            10000,  # Stem Cell Supporting Factors
            1,      # DAMPs (Tissue Damage)
            500,      # Activated leukocytes
            500,   # Stable leukocytes
            500,       # Suppressor leukocytes
            1,      # MDSC
            5000,  # MF

        ]

        output_keys = ['HQ', 'HM', 'N', 'P', 'A', 'SCSF', 'K', 'Q', 'S', 'U', 'MDSC', 'MF']

        init_dict = dict(zip(output_keys, init_state))      # this dict will be saved as a .csv if save_init_states = True

    if save_init_states:        # save initial values to a .csv preset   
        init_fname = path / f'init_val_preset_{run_number}.csv'

        with open(init_fname, "w", newline='') as file:
            writer = csv.writer(file)
            
            for key, value in init_dict.items():
                writer.writerow([key, value])
        
        print("Initial values successfully saved to .csv")

else:               # set initial values by reading from .csv

    # Warning: if .csv key-values do not match selected model format, error in sim will occur

    with open(read_from_init_loc, mode="r", encoding="utf-8") as file:

        reader = csv.reader(file)

        init_from_csv = {rows[0]: float(rows[1]) for rows in reader}

    if model == 2:
        
        init_state = [
            init_from_csv['HQ'],
            init_from_csv['HM'],
            init_from_csv['N'],
            init_from_csv['P'],
            init_from_csv['A'],
            init_from_csv['SCSF'],
            init_from_csv['K'],
            init_from_csv['Q'],
            init_from_csv['S'],
            init_from_csv['U']
        ]
    
    elif model == 3:
    
        init_state = [
            init_from_csv['HQ'],
            init_from_csv['HM'],
            init_from_csv['N'],
            init_from_csv['P'],
            init_from_csv['A'],
            init_from_csv['SCSF'],
            init_from_csv['K'],
            init_from_csv['Q'],
            init_from_csv['S'],
            init_from_csv['U'],
            init_from_csv['MDSC'],
            init_from_csv['MF']
        ]

    print("Initial values successfully read from .csv")


# parameters here
if not read_from_param:         # set model parameters manually
    
    if model == 2:
        parameters = {

            'k_H' : 3,
            'dH' : 0.05,
            'theta_N' : 100_000,       # 2000
            'theta_K' : 50_000,       # 5000
            'tau_Q' : 1,
            'tau_U' : 1,
            'd_SCSF' : 0.3,
            'd_S' : 0.7,
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
            'S_AS' : 3,
            'S_SCSF' : 10000,
            'S_KD' : 1,
            'k_sn' : 3,
            'k_nq' : 10,
            'k_ns' : 0.5,
            'R_KU' : 10,
            'I_crit' : 0.8,
            'K_crit' : 30_000,
            'k' : 1,
            'A_crit' : 3,
            'psi' : 1/10

        }

    elif model == 3:
        parameters = {

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

    if save_parameters:

        param_fname = path / f'parameter_preset_{run_number}.csv'

        with open(param_fname, "w", newline='') as file:
            writer = csv.writer(file)
            
            for key, value in parameters.items():
                writer.writerow([key, value])
        
        print("Model parameters successfully saved to .csv")

else:           # set model parameters by reading from .csv

    # Warning: if .csv key-values do not match selected model format, error in sim will occur
    with open(read_from_param_loc, mode="r", encoding="utf-8") as file:

        reader = csv.reader(file)

        parameters = {rows[0]: float(rows[1]) for rows in reader}

    print("Model parameters successfully read from provided .csv")
    


# ----- 0. Generating arrays before running solver -------------
if model == 2:
    num_outputs = 11
    ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

elif model == 3:
    num_outputs = 13
    ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

ext_stimuli = np.zeros((runs, num_outputs-1, len(timesteps)))

for i in range(runs):       # add stimuli here

    if delayed_infection == False:                          # one-time pathogen input
        ext_stimuli[i, 2, int(t_infection_start/delta_t)] = path_size_default + path_increment*i
        ext_stimuli[i, 2, int(nosocomial_start/delta_t)] = nosocomial_size        # optional; nosocomial infection

    else:                                                   # delayed input
        for j in np.arange(int(t_infection_start/delta_t), int(t_infection_end/delta_t)):

            ext_stimuli[i, 2, j] = (path_size_default + path_increment*i) / len(np.arange(int(t_infection_start/delta_t), int(t_infection_end/delta_t)))

        for j in np.arange(int(nosocomial_start/delta_t), int(nosocomial_end/delta_t)):

            ext_stimuli[i, 2, j] = (nosocomial_size) / len(np.arange(int(nosocomial_start/delta_t), int(nosocomial_end/delta_t)))

if HSPC_boost:

    for i in range(runs):

        ext_stimuli[i, 0, int(HSPC_boost_time/delta_t)] = HSPC_default + HSPC_increment*i

if SCSF_boost:

    for i in range(runs):

        ext_stimuli[i, 5, int(SCSF_boost_time/delta_t)] = SCSF_default + SCSF_increment*i

'''if bDerivatives:
    derivatives = np.zeros((runs, num_outputs-1, len(timesteps)))'''

outputs = np.zeros((runs, num_outputs, len(timesteps)))



# ------- 1. Running ODE solver (can easily swap this out for something like scpipy's solve_ivp() ) -------

start = time.time()
for i in range(runs):

    if model == 2:
        data = PL.lin_sim(PL.model_2_derivatives, parameters, init_state, t_final, delta_t, ext_stimuli[i], ext_stim_m)
        outputs[i, :, :] = data[0]
        print(f"Run {i+1} output successfully computed")

        '''if bDerivatives:
            derivatives[i] = data[1]
            print(f"Run {i+1} derivatives successfully loaded")'''
    
    elif model == 3:
        data = PL.lin_sim(beta_model_3, parameters, init_state, t_final, delta_t, ext_stimuli[i], ext_stim_m)
        outputs[i, :, :] = data[0]
        print(f"Run {i+1} output successfully computed")

        '''if bDerivatives:
            derivatives[i] = data[1]
            print(f"Run {i+1} derivatives successfully loaded")'''

end = time.time()
print(f'Execution successful. Time elapsed: {end-start}s')

# ---- 2. Save output to .csv's and plot the simulations -------
path = Path.cwd() / 'Runs'

if not Path.exists(path):
    Path.mkdir(path)

if model == 2:
    path = path / 'Model 2'
elif model == 3:
    path = path / 'Model 3'

if not Path.exists(path):
    Path.mkdir(path)

path = path / f'Exp_{run_number}'

if not Path.exists(path):
    Path.mkdir(path)

if output_to_csv:
    print("Saving outputs to .csv files now.")
    start = time.time()

    for i in range(runs):

        if model == 2:
            df = pd.DataFrame(outputs[i].T, columns=['HQ', 'HM', 'N', 'P', 'A', 'SCSF', 'K', 'Q', 'S', 'U', 'I'], index=list(map(str, timesteps)))
            df.to_csv(path / f'SIM_{run_number}_{i}_output.csv')

        elif model == 3:
            df = pd.DataFrame(outputs[i].T, columns=['HQ', 'HM', 'N', 'P', 'A', 'SCSF', 'K', 'Q', 'S', 'U', 'MDSC', 'MF', 'I'], index=list(map(str, timesteps)))
            df.to_csv(path / f'SIM_{run_number}_{i}_output.csv')

    end = time.time()
    print(f"Outputs successfully saved to .csv's. Time elapsed: {end-start}s")

# instantiate fig, axs objects here
if model == 2:
    fig1, axs1 = plt.subplots(3, 1)
    fig2, axs2 = plt.subplots(3, 1)
    fig3, axs3 = plt.subplots(3, 1)
    fig4, axs4 = plt.subplots(2, 1)
    figs = [fig1, fig2, fig3, fig4]
    axes = [axs1, axs2, axs3, axs4]

elif model == 3:
    fig1, axs1 = plt.subplots(3, 1)
    fig2, axs2 = plt.subplots(3, 1)
    fig3, axs3 = plt.subplots(3, 1)
    fig4, axs4 = plt.subplots(3, 1)
    fig5, axs5 = plt.subplots(1, 1)
    figs = [fig1, fig2, fig3, fig4, fig5]
    axes = [axs1, axs2, axs3, axs4, axs5]

# graph titles here
if model == 2:
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
        '$I(t)$'
    ]
elif model == 3:
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

# axs[i%4, i//4]

# plotting code here
# the if-else statements are used to ensure only one set of labels is produced per figure
if model == 2:
    for i in range(num_outputs):

        for j in range(runs):   # sim runs are split into 4 separate figures to make it more readable

            if i < 3:
                if i == 0:
                    axs1[i%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')      # labels only created once per figure
                else:
                    axs1[i%3].plot(timesteps, outputs[j, i])

                axs1[i%3].title.set_text(titles[i])
            
            elif i < 6:
                if i == 3:
                    axs2[(i-3)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')
                else:
                    axs2[(i-3)%3].plot(timesteps, outputs[j, i])

                axs2[(i-3)%3].title.set_text(titles[i])
            
            elif i < 9:
                if i == 6:
                    axs3[(i-6)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')
                else:
                    axs3[(i-6)%3].plot(timesteps, outputs[j, i])
                                       
                axs3[(i-6)%3].title.set_text(titles[i])
            
            else:
                if i ==9:
                    axs4[(i-9)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')
                else:
                    axs4[(i-9)%3].plot(timesteps, outputs[j, i])
                                       
                axs4[(i-9)%3].title.set_text(titles[i])

elif model == 3:
    for i in range(num_outputs):

        for j in range(runs):   # sim runs are split into 5 separate figures to make it more readable

            if i < 3:
                if i == 0:
                    axs1[i%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')
                else:
                    axs1[i%3].plot(timesteps, outputs[j, i])

                axs1[i%3].title.set_text(titles[i])
            
            elif i < 6:
                if i == 3:
                    axs2[(i-3)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')
                else:
                    axs2[(i-3)%3].plot(timesteps, outputs[j, i])
                
                axs2[(i-3)%3].title.set_text(titles[i])
            
            elif i < 9:
                if i == 6:
                    axs3[(i-6)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')
                else:
                    axs3[(i-6)%3].plot(timesteps, outputs[j, i])

                axs3[(i-6)%3].title.set_text(titles[i])
            
            elif i < 12:
                if i == 9:
                    axs4[(i-9)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')
                else:
                    axs4[(i-9)%3].plot(timesteps, outputs[j, i])

                axs4[(i-9)%3].title.set_text(titles[i])
            
            else:
                axs5.plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(t_infection_start/delta_t)]}; {delayed_infection}')

                axs5.title.set_text(titles[i])
                

# ----- misc options for tuning graphs as need arises ------
# create the y=0 lines here
if model == 2:
    axs1[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs1[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs1[2].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')

    axs2[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs2[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs2[2].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')

    axs3[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs3[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs3[2].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')

    axs4[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs4[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')

elif model == 3:
    axs1[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs1[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs1[2].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')

    axs2[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs2[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs2[2].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')

    axs3[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs3[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs3[2].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')

    axs4[0].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs4[1].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    axs4[2].plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')
    
    axs5.plot(timesteps, np.zeros(len(timesteps)), color=(0, 0, 0, 0.5), linestyle='--')


for fig in figs:
    fig.subplots_adjust(left=0.15, hspace=0.4)
    fig.set_size_inches(14, 9)
    fig.legend(loc="upper left", fontsize=6)

axs1[2].set_ylim((0, 100_000))        # optional; set y-limit for pathogen graph (useful if N(t) grows to carrying capacity)
#axs3[0].set_ylim((0, 50_000))        # optional; set y-limit for pathogen graph (useful if K(t) grows large)

# ----------------------------------------------------------
if model == 2:
    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()

    fig1.savefig(path / f'sim_{run_number}_Hq_Hp_N.png', dpi=300)
    fig2.savefig(path / f'sim_{run_number}_P_A_SCSF.png', dpi=300)
    fig3.savefig(path / f'sim_{run_number}_K_Q_S.png', dpi=300)
    fig4.savefig(path / f'sim_{run_number}_U_I.png', dpi=300)

elif model == 3:
    '''fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()'''

    fig1.savefig(path / f'sim_{run_number}_Hq_Hp_N.png', dpi=300)
    fig2.savefig(path / f'sim_{run_number}_P_A_SCSF.png', dpi=300)
    fig3.savefig(path / f'sim_{run_number}_K_Q_S.png', dpi=300)
    fig4.savefig(path / f'sim_{run_number}_U_MDSC_MF.png', dpi=300)
    fig5.savefig(path / f'sim_{run_number}_I.png', dpi=300)

plt.show()

