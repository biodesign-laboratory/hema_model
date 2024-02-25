import numpy as np
import matplotlib.pyplot as plt

from function import simulate

include_all = True
runs = 10

g_N = 0.2
K_PS = 0.7
K_AQ = 0.7
K_AS = 0.7
K_PU = 0.7
k_nq = 0.85
k_ns = 0.2
k_tn = 0.33
w = 0.0005
P_crit = 2500
A_crit = 1000
t_half_leukocytes = 7
t_half_P = 4.1375
t_half_A = 7
t_double_P = 15
gamma = 0.0000001
S_a = 1
S_n = 1
N_inf = 20000000
S_PQ = 0.33
S_PH = 0.01
S_PS = 0.02
S_AS = 0.04
S_AH = 0.01
S_AU = 0.33
theta_ps = 10_000_000
theta_ar = 10_000_000
theta_AS = 10_000_000
theta_UP = 10_000_000
Immune_start = 500
Active_start = 500
Immune_crit = 2500
Active_crit = 2500
y = 0.0008
d_s = 1/70
d_p = 1/4
d_a = 1/4
d_q = 1/4
d_u = 0.05

RATES = [g_N, K_PS, K_AQ, K_AS, K_PU, k_nq, k_ns, k_tn, w, P_crit, A_crit, S_a, S_n, N_inf, S_PQ, S_PH, S_PS, S_AS, S_AH, S_AU, 
         theta_ps, theta_ar, theta_AS, theta_UP, Immune_start, Active_start, Immune_crit, Active_crit, y, d_s, 
         d_p, d_a, d_q, d_u]

TFinal = 200
timestep_size = 0.01
t = np.arange(0, TFinal, timestep_size)

path_repeat_ts = [50,100,150]
path_repeat_size = [3000,3000,3000]

H_0 = 1000
N_0 = 0
P_0 = 600
A_0 = 1160
T_0 = 7000
a_0 = 0.00
b_0 = 1
e_0 = 0.00
E_star_0 = 1

n_steady = (-1*g_N + np.sqrt((g_N**2)-4*(-k_ns*7000)*(g_N/N_inf)))/(-2*g_N/N_inf)
            
'''
n_steady represents one of two possible theoretical steady state values with the assumption
that Q=0
Solved for by setting derivative equal to 0 and using the quadradtic formula
'''
N_vals = np.zeros(runs)

plot_label_nums = np.zeros(runs)

for i in range(runs):
    
    if i >= 1:
        N_vals[i] = N_0 + i*100 + 4000
        #plot_label_nums[i] = N_0 + i*100 + 4000
    else:
        N_vals[i] = N_0
        #plot_label_nums[i] = int(N_0)

init_vals = np.eye(runs, 9)

for i in range(runs):
    init_vals[i] = np.array([H_0, N_vals[i], T_0, a_0, b_0, e_0, P_0, A_0, E_star_0])



outputs = np.zeros((runs, 9, int(TFinal//timestep_size+1)))

for i in range(runs):
    outputs[i]=simulate(init_vals[i], RATES, timestep_size, TFinal, path_repeat_ts, path_repeat_size)



plot_titles = [ 'HSPC\'s' , 'External Force', 'Pro-inflammatory cytokines', 'Anti-inflammatory cytokines',
                 'Stable WBCs', 'Activated WBCs', 'Immunosuppressive WBCs', 'Total WBCs', 'Total Inflammation'
                ]


# ---------------- plots options -------------------
num_graphs = 9
set_path_ylim = True        # useful if any runs lead to out-of-control pathogen growth
ylim = 5000

# -- Base plots for each output from 0 <= t <= TFinal --
for i in range(num_graphs):
    
    plt.figure(i)
    for output in outputs:
        
        plt.plot(t, output[i], label="N_0: " + str(output[1][0]))
    
    plt.title(plot_titles[i])
    plt.legend(fontsize="10", loc="center left", bbox_to_anchor=(1, 0.5))

    if set_path_ylim and (i==1 or i==8):
        plt.ylim((0, ylim))

#plt.show()


# -- Plots focused on pathogen insults --
generate_zoomed_plots = True
radius_of_zoom = 20

if generate_zoomed_plots:

    for timestep in path_repeat_ts:
        for i in range(num_graphs*len(path_repeat_ts)):
        
            plt.figure(i)
            for output in outputs:

                    plt.plot(np.arange(timestep-radius_of_zoom, timestep+radius_of_zoom, step=timestep_size), output[i%len(path_repeat_ts)][int((timestep-radius_of_zoom)/timestep_size):int((timestep+radius_of_zoom)/timestep_size)], label="N_0: " + str(output[1][0]))
            
            plt.title(plot_titles[i%len(path_repeat_ts)])
            plt.legend(fontsize="10", loc="center left", bbox_to_anchor=(1, 0.5))

            if set_path_ylim and (i%len(path_repeat_ts)==1 or i%len(path_repeat_ts)==8):
                plt.ylim((0, ylim))

        plt.show()

