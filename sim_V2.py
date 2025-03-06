import matplotlib.pyplot as plt
import numpy as np
import project_library as PL
import time
from M3_beta import beta_model_3
import pandas as pd
from pathlib import Path

import matplotlib as mpl
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{siunitx}')
plt.rc('pgf.texsystem.pdflatex')

run_number = 'debug'                # used in file names, doesn't have to be a number
model = 3                           # 2 - model 2 (previous, no MDSCs), 3 - model 3 (MDSCs)

runs = 10
delta_t = 0.01
t_final = 672                       # 672 hours = 4 weeks (timescale is arbitrary but I still use it to gauge model behavior)
# num_outputs = 11
bDerivatives = True
graph_derivatives = True
delayed_infection = False            # use this to make pathogen input spread out over time period [t_infection_start, t_infection_end)
t_infection_start = 100             
t_infection_end = 115               # only used if delayed_infection = True
path_size_default = 10000
nosocomial_size = 20000
path_increment = 1000
output_to_csv = True

timesteps = np.arange(stop=t_final, step=delta_t)
print(len(timesteps))

# initial states here
if model == 2:
    init_state = [

        20000,  # Quiescent HSPCs
        0,      # Proliferating HSPCs
        0,      # PAMPs (Pathogens)
        10000,      # Pro-inflammatory Cytokines
        10000,      # Anti-inflammatory Cytokines
        10000,  # Stem Cell Supporting Factors
        0,      # DAMPs (Tissue Damage)
        0,      # Activated leukocytes
        1,   # Stable leukocytes
        0       # Suppressor leukocytes

    ]
elif model == 3:
    init_state = [

        20000,  # Quiescent HSPCs
        1,      # Proliferating HSPCs
        1,      # PAMPs (Pathogens)
        10000,      # Pro-inflammatory Cytokines
        10000,      # Anti-inflammatory Cytokines
        10000,  # Stem Cell Supporting Factors
        1,      # DAMPs (Tissue Damage)
        1,      # Activated leukocytes
        1,   # Stable leukocytes
        1,       # Suppressor leukocytes
        1,      # MDSC
        10000,  # MF

    ]

# parameters here
if model == 2:
    parameters = {

        'k_H' : 3,
        'dH' : 0.05,
        'theta_N' : 2000,
        'theta_K' : 5000,
        'tau_Q' : 1,
        'tau_U' : 1,
        'd_SCSF' : 0.3,
        'd_S' : 0.7,
        'd_Q' : 0.9,
        'd_U' : 0.5,
        'd_P' : 0.95,
        'd_A' : 0.95,
        'g_N' : 0.10,
        'N_oo' : 2 * 10**7,
        'N_half' : 2500,
        'S_PH' : 5,
        'S_PS' : 5,
        'S_PQ' : 20,
        'S_AU' : 15,
        'S_AH' : 5,
        'S_AS' : 5,
        'S_SCSF' : 20000,
        'S_KD' : 1,
        'k_sn' : 3,
        'k_nq' : 10,
        'k_ns' : 0.5,
        'R_KU' : 10,
        'I_crit' : 0.8,
        'K_crit' : 10000,
        'k' : 1,
        'A_crit' : 3,
        'psi' : 1/10

    }
elif model == 3:
    parameters = {

        'k_H' : 3,
        'dH' : 0.05,
        'theta_N' : 2000,
        'theta_K' : 5000,
        'tau_Q' : 1,
        'tau_U' : 1,
        'd_SCSF' : 0.3,
        'd_S' : 0.7,
        'd_Q' : 0.7,
        'd_U' : 0.5,
        'd_P' : 0.95,
        'd_A' : 0.95,
        'g_N' : 0.10,
        'N_oo' : 2 * 10**7,
        'N_half' : 2500,
        'S_PH' : 2,
        'S_PS' : 4,
        'S_PQ' : 10,
        'S_AU' : 15,
        'S_AH' : 0,
        'S_AS' : 3,
        'S_AM' : 7,
        'S_SCSF' : 10000,
        'S_KD' : 1,
        'k_sn' : 3,
        'k_nq' : 10,
        'k_ns' : 0.5,
        'R_KU' : 10,
        'I_crit' : 0.8,
        'K_crit' : 100_000,
        'k' : 1,
        'psi' : 1,
        'd_M' : 9/10,
        'd_MF' : 0.3,
        'S_KMD' : 1/5,
        'S_KQ' : 1/3,
        'C_QM' : 1,
        'C_MDM' : 1/2,
        'C_UM' : 1/3,
        'S_MF' : 10000

    }

# ----- 0. Generating arrays before running solver -------------
if model == 2:
    num_outputs = 11
    ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

elif model == 3:
    num_outputs = 13
    ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

ext_stimuli = np.zeros((runs, num_outputs-1, len(timesteps)))

for i in range(runs):       # add stimuli here

    if delayed_infection == False:                          # big, one-time pathogen input
        ext_stimuli[i, 2, int(t_infection_start/delta_t)] = path_size_default + path_increment*i
        ext_stimuli[i, 2, int(300/delta_t)] = path_size_default        # optional; nosocomial infection

    else:                                                   # delayed input
        for j in np.arange(int(t_infection_start/delta_t), int(t_infection_end/delta_t)):

            ext_stimuli[i, 2, j] = (path_size_default + path_increment*i) / len(np.arange(int(t_infection_start/delta_t), int(t_infection_end/delta_t)))

if bDerivatives:
    derivatives = np.zeros((runs, num_outputs-1, len(timesteps)))

outputs = np.zeros((runs, num_outputs, len(timesteps)))



# ------- 1. Running ODE solver -------

start = time.time()
for i in range(runs):

    if model == 2:
        data = PL.lin_sim(PL.model_2_derivatives, parameters, init_state, t_final, delta_t, ext_stimuli[i], ext_stim_m, return_derivatives=bDerivatives)
        outputs[i, :, :] = data[0]
        print(f"Run {i+1} output successfully computed")

        if bDerivatives:
            derivatives[i] = data[1]
            print(f"Run {i+1} derivatives successfully loaded")
    
    elif model == 3:
        data = PL.lin_sim(beta_model_3, parameters, init_state, t_final, delta_t, ext_stimuli[i], ext_stim_m, return_derivatives=bDerivatives)
        outputs[i, :, :] = data[0]
        print(f"Run {i+1} output successfully computed")

        if bDerivatives:
            derivatives[i] = data[1]
            print(f"Run {i+1} derivatives successfully loaded")

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

elif model == 3:
    fig1, axs1 = plt.subplots(3, 1,figsize=(8,6))
    fig2, axs2 = plt.subplots(3, 1,figsize=(8,6))
    fig3, axs3 = plt.subplots(3, 1,figsize=(8,6))
    fig4, axs4 = plt.subplots(3, 1,figsize=(8,6))
    fig5, axs5 = plt.subplots(1, 1,figsize=(8,6))

# graph titles here
if model == 2:
    titles = [
        '$H_Q(t)$',
        '$H_P(t)$',
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
        '$H_P(t)$',
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
if model == 2:
    for i in range(num_outputs):

        for j in range(runs):   # sim runs are split into 4 separate figures to make it more readable

            if i < 3:
                axs1[i%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
                axs1[i%3].title.set_text(titles[i])
                axs1[i%3].legend()
            
            elif i < 6:
                axs2[(i-3)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
                axs2[(i-3)%3].title.set_text(titles[i])
                axs2[i%3].legend()
            
            elif i < 9:
                axs3[(i-6)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
                axs3[(i-6)%3].plot(timesteps, np.zeros(len(outputs[j, i])))
                axs3[(i-6)%3].title.set_text(titles[i])
                axs3[i%3].legend()
            
            else:
                axs4[(i-9)%3].plot(timesteps, outputs[j, i], label=f'{ext_stimuli[j, 2, int(100/delta_t)]}')
                axs4[(i-9)%3].plot(timesteps, np.zeros(len(outputs[j, i])))
                axs4[(i-9)%3].title.set_text(titles[i])
                axs4[i%3].legend()

elif model == 3:
    for i in range(num_outputs):

        for j in range(runs):   # sim runs are split into 4 separate figures to make it more readable

            if i < 3:
                axs1[i%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
                axs1[i%3].title.set_text(titles[i])
                axs1[i%3].legend()
            
            elif i < 6:
                axs2[(i-3)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
                axs2[(i-3)%3].title.set_text(titles[i])
                axs2[i%3].legend()
            
            elif i < 9:
                axs3[(i-6)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
                axs3[(i-6)%3].plot(timesteps, np.zeros(len(outputs[j, i])))
                axs3[(i-6)%3].title.set_text(titles[i])
                axs3[i%3].legend()
            
            elif i < 12:
                axs4[(i-9)%3].plot(timesteps, outputs[j, i], label=f'{ext_stimuli[j, 2, int(100/delta_t)]}')
                axs4[(i-9)%3].plot(timesteps, np.zeros(len(outputs[j, i])))
                axs4[(i-9)%3].title.set_text(titles[i])
                axs4[i%3].legend()
            
            else:
                axs5.plot(timesteps, outputs[j, i], label=f'{ext_stimuli[j, 2, int(100/delta_t)]}')
                axs5.plot(timesteps, np.zeros(len(outputs[j, i])))
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

axs1[2].set_ylim((0, 150_000))        # optional; set y-limit for pathogen graph (useful if N(t) grows to carrying capacity)
#axs2[2].set_ylim((0, 150_000))        # optional; set y-limit for pathogen graph (useful if SCSF(t) grows to carrying capacity)

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

if model == 3:
    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()

    fig1.savefig(path / f'sim_{run_number}_Hq_Hp_N.png', dpi=300)
    fig2.savefig(path / f'sim_{run_number}_P_A_SCSF.png', dpi=300)
    fig3.savefig(path / f'sim_{run_number}_K_Q_S.png', dpi=300)
    fig4.savefig(path / f'sim_{run_number}_U_MDSC_MF.png', dpi=300)
    fig5.savefig(path / f'sim_{run_number}_I.png', dpi=300)

#plt.show()

# ------------- 3. Saving derivatives to .csv files -------------

if bDerivatives:

    print("Saving derivatives to .csv files now.")
    start = time.time()

    for i in range(runs):

        df = pd.DataFrame(derivatives[i].T, columns=['dHQ', 'dHM', 'dN', 'dP', 'dA', 'dSCSF', 'dK', 'dQ', 'dS', 'dU'], index=list(map(str, timesteps)))
        df.to_csv(path / f'SIM_{run_number}_{i}_derivatives.csv')
    
    end = time.time()
    print(f'Derivative data successfully saved. Time elapsed: {end-start}')
