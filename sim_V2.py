import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import project_library as PL
import time
# import M2_beta
import pandas as pd
from pathlib import Path
from matplotlib.collections import LineCollection
import global_variables as gv


run_number = 'debug'      # used in file names

runs = 10
experiments = 10
delta_t = 0.01
t_final = 672       # 672 hours = 4 weeks
num_outputs = 11
bDerivatives = False
bDebug = False
graph_derivatives = False

timesteps = np.arange(stop=t_final, step=delta_t)
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

ext_stimuli = np.zeros((runs, num_outputs, len(timesteps)))

if bDerivatives:
    derivatives = np.zeros((runs, 10, len(timesteps)))

outputs = np.zeros((runs, num_outputs, len(timesteps)))

ext_stim_m = ['ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD', 'ADD']

key_names = [
    'k_H', 
    'dH', 
    'theta_N',
    'theta_K',
    'tau_Q',
    'tau_U',
    'd_SCSF',
    'd_S',
    'd_Q',
    'd_U',
    'd_P',
    'd_A',
    'g_N',
    'N_oo',
    'N_half',
    'S_PH',
    'S_PS',
    'S_PQ',
    'S_AU',
    'S_AH',
    'S_AS',
    'S_SCSF',
    'S_KD',
    'k_sn',
    'k_nq',
    'k_ns',
    'R_KU',
    'I_crit',
    'K_crit',
    'k',
    'A_crit',
    'psi']

np.random.seed(0)

start = time.time()
sol_total = [[0] * runs] * experiments

fig1, axs1 = plt.subplots(3, 1)
fig2, axs2 = plt.subplots(3, 1)
fig3, axs3 = plt.subplots(3, 1)
fig4, axs4 = plt.subplots(2, 1)

for key in key_names:
    parameter_range = np.linspace(gv.parameters_global[key] / 4, (gv.parameters_global[key] * 7) / 4, experiments)
    for h in range(experiments):
        gv.parameters_global[key] = parameter_range[h]
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
            sol = solve_ivp(PL.model_2_derivatives, [0, t_final], init_state)
            print(f"Run {i} output successfully computed")
            sol_total[h][i] = sol

        # ---- plot the simulations and save output to png's -------
        path = Path.cwd() / 'Experiments_1'
        
        if not Path.exists(path):
            Path.mkdir(path)

        path = path / f'Runs_{key}'

        if not Path.exists(path):
            Path.mkdir(path)

        path = path / f'Exp_{h}'

        if not Path.exists(path):
            Path.mkdir(path)

        # ---- plotting the simulations -------

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
                sol = sol_total[h][i]
                if j < 3:
                    # axs1[j%3].set_ylim((0, 100000))
                    axs1[j%3].plot(sol.t, sol.y[j - 1])
                    axs1[j%3].autoscale(enable=True, axis='y')
                    axs1[j%3].set_xlim((0, t_final))
                    axs1[j%3].title.set_text(titles[j - 1])
                    # axs1[j%3].plot(sol.t, sol.y[j - 1])
                    # axs1[j%3].add_collection(lc[j])
                    
                elif j < 6:
                    # axs2[(j-3)%3].set_ylim((0, 100000))
                    axs2[(j-3)%3].plot(sol.t, sol.y[j - 1])
                    axs2[(j-3)%3].autoscale(enable=True, axis='y')
                    axs2[(j-3)%3].set_xlim((0, t_final))
                    axs2[(j-3)%3].title.set_text(titles[j - 1])
                    # axs2[(j-3)%3].plot(sol.t, sol.y[j - 1])
                    # axs2[(j-3)%3].add_collection(lc[j])
                
                elif j < 9:
                    # axs3[(j-6)%3].set_ylim((0, 100000))
                    axs3[(j-6)%3].plot(sol.t, sol.y[j - 1])
                    axs3[(j-6)%3].autoscale(enable=True, axis='y')
                    axs3[(j-6)%3].set_xlim((0, t_final))
                    axs3[(j-6)%3].title.set_text(titles[j - 1])
                    # axs3[(j-6)%3].plot(sol.t, sol.y[j - 1])
                    # axs3[(j-6)%3].add_collection(lc[j])
                
                else:
                    # axs4[(j-9)%3].set_ylim((0, 100000))
                    axs4[(j-9)%3].plot(sol.t, sol.y[j - 1])
                    axs4[(j-9)%3].autoscale(enable=True, axis='y')
                    axs4[(j-9)%3].set_xlim((0, t_final))
                    axs4[(j-9)%3].title.set_text(titles[j - 1])
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

        fig1.savefig(path / f'exp_{h}_Hq_Hp_N.png', dpi=300)
        fig2.savefig(path / f'exp_{h}_P_A_SCSF.png', dpi=300)
        fig3.savefig(path / f'exp_{h}_K_Q_S.png', dpi=300)
        fig4.savefig(path / f'exp_{h}_U_I.png', dpi=300)
        print(f"Experiment {h} successfully saved")

        for i in range(3):
            axs1[i].clear()
            axs2[i].clear()
            axs3[i].clear()

        for i in range(2):
            axs4[i].clear()

    end = time.time()
    print(f'Execution successful. Time elapsed: {end-start}s')
