import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import project_library as PL
import time
# import M2_beta
import pandas as pd
from pathlib import Path
from matplotlib.collections import LineCollection


run_number = 'debug'      # used in file names

runs = 5000
delta_t = 0.01
t_final = 672       # 672 hours = 4 weeks
num_outputs = 11
bDerivatives = False
bDebug = False
graph_derivatives = False

timesteps = np.arange(stop=t_final, step=delta_t)
print(len(timesteps))
init_state = [

    11100,  # Quiescent HSPCs
    1000,      # Proliferating HSPCs
    1,      # PAMPs (Pathogens)
    1000,      # Pro-inflammatory Cytokines
    1000,      # Anti-inflammatory Cytokines
    32000,  # Stem Cell Supporting Factors
    0,      # DAMPs (Tissue Damage)
    1,      # Activated leukocytes
    10,   # Stable leukocytes
    0.5       # Suppressor leukocytes

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
    'S_PQ' : 5,
    'S_AU' : 7,
    'S_AH' : 3,
    'S_AS' : 1,
    'S_SCSF' : 10000,
    'S_KD' : 1,
    'k_sn' : 3,
    'k_nq' : 10,
    'k_ns' : 0.5,
    'R_KU' : 15,
    'I_crit' : 0.8,
    'K_crit' : 10000,
    'k' : 1,
    'A_crit' : 3,
    'psi' : 1/10

}

ext_stimuli = np.zeros((runs, num_outputs, len(timesteps)))

if bDerivatives:
    derivatives = np.zeros((runs, 10, len(timesteps)))

outputs = np.zeros((runs, num_outputs, len(timesteps)))

ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

start = time.time()
sol_total = [0] * runs
for i in range(runs):
    init_state = [

    40000 * np.random.rand(),  # Quiescent HSPCs
    2000 * np.random.rand(),      # Proliferating HSPCs
    100000 * np.random.rand(),      # PAMPs (Pathogens)
    20000 * np.random.rand(),      # Pro-inflammatory Cytokines
    20000 * np.random.rand(),      # Anti-inflammatory Cytokines
    20000 * np.random.rand(),  # Stem Cell Supporting Factors
    100000 * np.random.rand(),      # DAMPs (Tissue Damage)
    2 * np.random.rand(),      # Activated leukocytes
    20 * np.random.rand(),   # Stable leukocytes
    1 * np.random.rand()       # Suppressor leukocytes

]
    sol = solve_ivp(PL.model_x_derivatives, [0, t_final], init_state)
    print(f"Run {i} output successfully computed")
    sol_total[i] = sol


# plt.plot(sol.t, sol.y[0])
# plt.xlabel('Time')
# plt.ylabel('y')
# plt.show()

end = time.time()
print(f'Execution successful. Time elapsed: {end-start}s')

# ---- plotting the simulations -------
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

lc = [0] * num_outputs
segments_total = [0] * runs

for j in range(num_outputs):
    for i in range(runs):
        # segments = outputs[i,j]
        # segments_total[i] = np.column_stack((timesteps, segments))

    # lc[j] = LineCollection(segments_total, linewidths=0.5, alpha=0.3, color='blue')
        sol = sol_total[i]
        if j < 3:
            # axs1[j%3].set_ylim((0, 100000))
            axs1[j%3].plot(sol.t, sol.y[j - 1])
            axs1[j%3].autoscale(enable=True, axis='y')
            axs1[j%3].set_xlim((0, t_final))
            axs1[j%3].title.set_text(titles[j])
            # axs1[j%3].plot(sol.t, sol.y[j - 1])
            # axs1[j%3].add_collection(lc[j])
            
        elif j < 6:
            # axs2[(j-3)%3].set_ylim((0, 100000))
            axs2[(j-3)%3].plot(sol.t, sol.y[j - 1])
            axs2[(j-3)%3].autoscale(enable=True, axis='y')
            axs2[(j-3)%3].set_xlim((0, t_final))
            axs2[(j-3)%3].title.set_text(titles[j])
            # axs2[(j-3)%3].plot(sol.t, sol.y[j - 1])
            # axs2[(j-3)%3].add_collection(lc[j])
        
        elif j < 9:
            # axs3[(j-6)%3].set_ylim((0, 100000))
            axs3[(j-6)%3].plot(sol.t, sol.y[j - 1])
            axs3[(j-6)%3].autoscale(enable=True, axis='y')
            axs3[(j-6)%3].set_xlim((0, t_final))
            axs3[(j-6)%3].title.set_text(titles[j])
            # axs3[(j-6)%3].plot(sol.t, sol.y[j - 1])
            # axs3[(j-6)%3].add_collection(lc[j])
        
        else:
            # axs4[(j-9)%3].set_ylim((0, 100000))
            axs4[(j-9)%3].plot(sol.t, sol.y[j - 1])
            axs4[(j-9)%3].autoscale(enable=True, axis='y')
            axs4[(j-9)%3].set_xlim((0, t_final))
            axs4[(j-9)%3].title.set_text(titles[j])
            # axs4[(j-9)%3].plot(sol.t, sol.y[j - 1])
            # axs4[(j-9)%3].add_collection(lc[j])


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


fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()

plt.show()
