import matplotlib.pyplot as plt
import numpy as np
import project_library as PL
import time
# import M2_beta
import pandas as pd
from pathlib import Path

run_number = 'debug'      # used in file names

runs = 5
delta_t = 0.01
t_final = 672       # 672 hours = 4 weeks
num_outputs = 11
bDerivatives = True
bDebug = True
graph_derivatives = True

timesteps = np.arange(stop=t_final, step=delta_t)
print(len(timesteps))
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
    'd_U' : 0.8,
    'd_P' : 0.95,
    'd_A' : 0.95,
    'g_N' : 0.10,
    'N_oo' : 2 * 10**7,
    'N_half' : 500,
    'S_PH' : 3,
    'S_PS' : 1,
    'S_PQ' : 10,
    'S_AU' : 7,
    'S_AH' : 3,
    'S_AS' : 1,
    'S_SCSF' : 10000,
    'S_KD' : 1,
    'k_sn' : 3,
    'k_nq' : 10,
    'k_ns' : 0.5,
    'R_KU' : 5,
    'I_crit' : 0.8,
    'K_crit' : 100000,
    'k' : 1,
    'A_crit' : 3,
    'psi' : 1/10

}

ext_stimuli = np.zeros((runs, num_outputs, len(timesteps)))

for i in range(runs):       # add stimuli here

    ext_stimuli[i, 2, int(100/delta_t)] = 20000 + 10000*i
    # ext_stimuli[i, 2, int(300/delta_t)] = 5000        # optional; nosocomial infection

if bDerivatives:
    derivatives = np.zeros((runs, 10, len(timesteps)))

if bDebug:
    debug_output = np.zeros((runs, 10, len(timesteps)))     # change 10 to however many terms you have returned by ODE_eq()

outputs = np.zeros((runs, num_outputs, len(timesteps)))

ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

start = time.time()
for i in range(runs):
    data = PL.lin_sim(PL.model_3_derivatives, parameters, init_state, t_final, delta_t, ext_stimuli[i], ext_stim_m, return_derivatives=bDerivatives)
    outputs[i, :, :] = data[0]
    print(f"Run {i} output successfully computed")

    if bDerivatives:
        derivatives[i] = data[1]
        print(f"Run {i} derivatives successfully loaded")
    if bDebug:
        debug_output = data[2]
        print(f"Run {i} debug data successfully loaded")

end = time.time()
print(f'Execution successful. Time elapsed: {end-start}s')

# ---- plotting the simulations -------
path = Path.cwd() / 'Runs'

if not Path.exists(path):
    Path.mkdir(path)

path = Path.cwd() / 'Runs' / f'Exp_{run_number}'

if not Path.exists(path):
    Path.mkdir(path)

'''path = Path.cwd() / 'Runs' / f'Exp_{run_number}' / f'sim_{run_number}'
if not Path.exists(path):
    Path.mkdir(path)'''


fig1, axs1 = plt.subplots(3, 1)
fig2, axs2 = plt.subplots(3, 1)
fig3, axs3 = plt.subplots(3, 1)
fig4, axs4 = plt.subplots(2, 1)

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

# axs[i%4, i//4]

for i in range(num_outputs):

    for j in range(runs):   # sim runs are split into 4 separate figures to make it more readable

        if i < 3:
            axs1[i%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
            axs1[i%3].title.set_text(titles[i])
            #axs1[i%3].legend()
        
        elif i < 6:
            axs2[(i-3)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
            axs2[(i-3)%3].title.set_text(titles[i])
            #axs2[i%3].legend()
        
        elif i < 9:
            axs3[(i-6)%3].plot(timesteps, outputs[j, i], label=f'N={ext_stimuli[j, 2, int(100/delta_t)]}')
            axs3[(i-6)%3].plot(timesteps, np.zeros(len(outputs[j, i])))
            axs3[(i-6)%3].title.set_text(titles[i])
            #axs3[i%3].legend()
        
        else:
            axs4[(i-9)%3].plot(timesteps, outputs[j, i], label=f'{ext_stimuli[j, 2, int(100/delta_t)]}')
            axs4[(i-9)%3].plot(timesteps, np.zeros(len(outputs[j, i])))
            axs4[(i-9)%3].title.set_text(titles[i])
            #axs4[i%3].legend()

# ----- misc options for tuning graphs as need arises ------

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

axs1[2].set_ylim((0, 70000))        # optional; set y-limit for pathogen graph if N(t) grows to carrying capacity

# ----------------------------------------------------------

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()

fig1.savefig(path / f'sim_{run_number}_Hq_Hp_N.png', dpi=300)
fig2.savefig(path / f'sim_{run_number}_P_A_SCSF.png', dpi=300)
fig3.savefig(path / f'sim_{run_number}_K_Q_S.png', dpi=300)
fig4.savefig(path / f'sim_{run_number}_U_I.png', dpi=300)

plt.show()

# ------------- saving derivatives and/or debug terms to .csv files -------------

if bDerivatives:

    print("Saving derivatives to .csv files now.")
    start = time.time()

    for i in range(runs):

        df = pd.DataFrame(derivatives[i].T, columns=['dHQ', 'dHM', 'dN', 'dP', 'dA', 'dSCSF', 'dK', 'dQ', 'dS', 'dU'], index=list(map(str, timesteps)))
        df.to_csv(path / f'SIM_{run_number}_{i}_derivatives.csv')
    
    end = time.time()
    print(f'Derivative data successfully saved. Time elapsed: {end-start}')

if bDebug:

    print("Saving debug data to .csv files now.")
    start = time.time()

    for i in range(runs):

        df = pd.DataFrame(derivatives[i].T, columns=['eta_Q', 'eta_M', 'D_I', 'beta', 'dH * H_M', 'D_Q', 'D_U', 'N kills S', 'gN*N*(1-N/Noo)', 'S + Q kills N'], index=list(map(str, timesteps)))
        df.to_csv(path / f'SIM_{run_number}_{i}_debug.csv')
    
    end = time.time()
    print(f'Debug data successfully saved. Time elapsed: {end-start}')

# ----------- graph derivatives -------------------
if graph_derivatives:

    fig, axs = plt.subplots(3, 4)
    titles = [
        '$dH_Q/dt$',
        '$dH_P/dt$',
        '$dN/dt$',
        '$dP/dt$',
        '$dA/dt$',
        '$dSCSF/dt$',
        '$dK/dt$',
        '$dQ/dt$',
        '$dS/dt$',
        '$dU/dt$', 
    ]

    for i in range(10):

        for j in range(runs):

            axs[i%3, i//3].plot(timesteps, derivatives[j, i])
            axs[i%3, i//3].title.set_text(titles[i])


    path = Path.cwd() / 'Runs' / f'Exp_{run_number}' / f'derivatives_{run_number}.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()

    # What needs to be done:
    #   - Split figures into separate graphs so png's are easier to read
