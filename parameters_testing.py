import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import project_library as PL
import time
from pathlib import Path
import global_variables as gv

runs = 500
experiments = 10
delta_t = 0.01
t_final = 672       # 672 hours = 4 weeks
delta_t = 0.01
num_outputs = 12

timesteps = np.arange(stop=t_final, step=delta_t)

key_names_2 = gv.key_names_2
key_names_3 = gv.key_names_3

np.random.seed(0)

init_state_3 = np.array([

            40000 * np.random.rand(runs),  # Quiescent HSPCs
            2000 * np.random.rand(runs),      # Proliferating HSPCs
            100000 * np.random.rand(runs),      # PAMPs (Pathogens)
            20000 * np.random.rand(runs),      # Pro-inflammatory Cytokines
            20000 * np.random.rand(runs),      # Anti-inflammatory Cytokines
            20000 * np.random.rand(runs),  # Stem Cell Supporting Factors
            100000 * np.random.rand(runs),      # DAMPs (Tissue Damage)
            2 * np.random.rand(runs),      # Activated leukocytes
            20 * np.random.rand(runs),   # Stable leukocytes
            1 * np.random.rand(runs),       # Suppressor leukocytes
            1 * np.random.rand(runs),      # MDSC
            5000 * np.random.rand(runs),  # MF
        ])

start = time.time()
sol_total = [[0] * runs] * experiments

fig1, axs1 = plt.subplots(3, 1)
fig2, axs2 = plt.subplots(3, 1)
fig3, axs3 = plt.subplots(3, 1)
fig4, axs4 = plt.subplots(3, 1)

# Script for running parameter testing
'''
for key in key_names_3:
    parameter_range = np.linspace(gv.parameters_3[key] / 8, (gv.parameters_3[key] * 15) / 8, experiments)
    print(parameter_range)
    start_key = time.time()
    if key == 'd_P':
            runs = 500
            experiments = 1
    for h in range(experiments):
        gv.parameters_3[key] = parameter_range[h]
        runs = 1
        if key == 'd_P':
            runs = 500
            gv.parameters_3[key] = 1.04236111

        for i in range(runs):
            sol = solve_ivp(PL.model_3_derivatives, [0, t_final], init_state_3[:, i])
            print(f"Run {i} output successfully computed")
            sol_total[h][i] = sol

        if key == 'd_P':
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
                '$MDSC(t)$',
                '$MF(t)$'
            ]

            for j in range(num_outputs):
                for i in range(runs):
                    sol = sol_total[h][i]
                    if j < 3:
                        axs1[j%3].plot(sol.t, sol.y[j])
                        axs1[j%3].autoscale(enable=True, axis='y')
                        axs1[j%3].set_xlim((0, t_final))
                        axs1[j%3].title.set_text(titles[j])
                        
                    elif j < 6:
                        axs2[(j-3)%3].plot(sol.t, sol.y[j])
                        axs2[(j-3)%3].autoscale(enable=True, axis='y')
                        axs2[(j-3)%3].set_xlim((0, t_final))
                        axs2[(j-3)%3].title.set_text(titles[j])
                    
                    
                    elif j < 9:
                        axs3[(j-6)%3].plot(sol.t, sol.y[j])
                        axs3[(j-6)%3].autoscale(enable=True, axis='y')
                        axs3[(j-6)%3].set_xlim((0, t_final))
                        axs3[(j-6)%3].title.set_text(titles[j])
                    
                    elif j < 12:
                        axs4[(j-9)%3].plot(sol.t, sol.y[j])
                        axs4[(j-9)%3].autoscale(enable=True, axis='y')
                        axs4[(j-9)%3].set_xlim((0, t_final))
                        axs4[(j-9)%3].title.set_text(titles[j])



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

            fig1.tight_layout()
            fig2.tight_layout()
            fig3.tight_layout()
            fig4.tight_layout()

            fig1.savefig(path / f'exp_{h}_Hq_Hp_N.png', dpi=300)
            fig2.savefig(path / f'exp_{h}_P_A_SCSF.png', dpi=300)
            fig3.savefig(path / f'exp_{h}_K_Q_S.png', dpi=300)
            fig4.savefig(path / f'exp_{h}_U_MDSC_MF.png', dpi=300)
            end_key = time.time()
            print(f"Experiment {h} successfully saved: {end_key - start_key}")

            for i in range(3):
                axs1[i].clear()
                axs2[i].clear()
                axs3[i].clear()
                axs4[i].clear()


end = time.time()
print(f'Execution successful. Time elapsed: {end-start}s')
'''

# Script for recreating oscillation runs
sol_total = [0] * 3
init_state_3 = [gv.ics_oscillating_1_1, gv.ics_oscillating_2, gv.ics_oscillating_3]
for i in range(3):
    sol = solve_ivp(PL.model_3_derivatives, [0, t_final], init_state_3[i], dense_output=True, rtol=1e-2, atol=1e-2, method='LSODA')
    sol_total = sol

    # ---- plot the simulations and save output to png's -------
    path = Path.cwd() / 'Experiments_2'

    if not Path.exists(path):
        Path.mkdir(path)

    path = path / f'Runs_Oscillation'

    if not Path.exists(path):
        Path.mkdir(path)

    path = path / f'Oscillation_{i}'

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
        '$MDSC(t)$',
        '$MF(t)$'
    ]

    for j in range(num_outputs):
        for i in range(3):
            if j < 3:
                axs1[j%3].plot(sol.t, sol.y[j])
                axs1[j%3].autoscale(enable=True, axis='y')
                axs1[j%3].set_xlim((0, t_final))
                axs1[j%3].title.set_text(titles[j])
                
            elif j < 6:
                axs2[(j-3)%3].plot(sol.t, sol.y[j])
                axs2[(j-3)%3].autoscale(enable=True, axis='y')
                axs2[(j-3)%3].set_xlim((0, t_final))
                axs2[(j-3)%3].title.set_text(titles[j])
            
            
            elif j < 9:
                axs3[(j-6)%3].plot(sol.t, sol.y[j])
                axs3[(j-6)%3].autoscale(enable=True, axis='y')
                axs3[(j-6)%3].set_xlim((0, t_final))
                axs3[(j-6)%3].title.set_text(titles[j])
            
            elif j < 12:
                axs4[(j-9)%3].plot(sol.t, sol.y[j])
                axs4[(j-9)%3].autoscale(enable=True, axis='y')
                axs4[(j-9)%3].set_xlim((0, t_final))
                axs4[(j-9)%3].title.set_text(titles[j])



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

    axs1[0].set_ylim((0, 2000)) 

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()

    fig1.savefig(path / f'exp_{0}_Hq_Hp_N.png', dpi=300)
    fig2.savefig(path / f'exp_{0}_P_A_SCSF.png', dpi=300)
    fig3.savefig(path / f'exp_{0}_K_Q_S.png', dpi=300)
    fig4.savefig(path / f'exp_{0}_U_MDSC_MF.png', dpi=300)


    for i in range(3):
        axs1[i].clear()
        axs2[i].clear()
        axs3[i].clear()
        axs4[i].clear()
print("Execution Successful")