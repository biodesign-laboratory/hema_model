import numpy as np
import matplotlib.pyplot as plt

from project_library import linear_sim
import time                                 # for measuring the execution time of linear_sim

include_all = True
runs = 1

g_N = 0.2
k_nq = 0.85
k_ns = 0.2
k_tn = 0.33
w = 0.0005
P_crit = 4000
t_half_leukocytes = 7
t_half_P = 4.1375
t_half_A = 7
t_double_P = 15
gamma = 0.0000001
S_a = 1
S_n = 1
N_inf = 20000000
S_PQ = 0.33
S_PH = 0.2
S_PS = 0.02
S_AS = 0.005
S_AH = 0.005
S_AU = 0.2
Immune_start = 1000
Active_start = 500
Immune_crit = 5000
Active_crit = 2500
y = 0.0008
d_s = 1/70
d_p = 1/4
d_a = 1/4
d_q = 1/10
d_u = 1/4

RATES = [g_N, k_nq, k_ns, k_tn, w, P_crit, S_PQ, S_PH, S_PS, S_AS, S_AH, S_AU, Immune_start, Active_start,
         Immune_crit, Active_crit, y, d_s, d_p, d_a, d_q, d_u]

TFinal = 200
timestep_size = 0.01
t = np.arange(0, TFinal, timestep_size)

path_repeat_ts = [100]      # which timesteps 
path_repeat_size = [5000]   # how strong / how many pathogens

path_repeat_arr = []
for n in range(runs):
    
    path_size = [5200 + 20*n]
    path_repeat_arr.append([path_repeat_ts, path_size])

H_0 = 1000
N_0 = 0
P_0 = 600
A_0 = 1160
T_0 = 2000
a_0 = 0.00
b_0 = 1
e_0 = 0.00
E_star_0 = 1

N_vals = np.zeros(runs)

plot_label_nums = np.zeros(runs)

for i in range(runs):
    
    if i >= 1:
        N_vals[i] = N_0 + i*200 + 2000

    else:
        N_vals[i] = N_0


init_vals = np.eye(runs, 9)

for i in range(runs):
    init_vals[i] = np.array([H_0, N_vals[i], T_0, a_0, b_0, e_0, P_0, A_0, E_star_0])



outputs = np.zeros((runs, 9, int(TFinal//timestep_size+1)))

start_time = time.time()
for i in range(runs):
    
    outputs[i]=linear_sim(init_vals[i], RATES, timestep_size, TFinal, path_repeat_arr[i][0], path_repeat_arr[i][1])

end_time = time.time()
run_time = end_time - start_time


plot_titles = [ 'HSPC\'s' , 'External Force', 'Pro-inflammatory cytokines', 'Anti-inflammatory cytokines',
                 'Stable WBCs', 'Activated WBCs', 'Immunosuppressive WBCs', 'Total WBCs', 'Total Inflammation'
                ]


# ---------------- plots: -------------------
num_graphs = 9
# -- Base plots for each output from 0 <= t <= TFinal --
for i in range(num_graphs):
    
    plt.figure(i)
    for output in outputs:
        
        plt.plot(t, output[i], label="path_size: " + str(int(output[1][int(100/timestep_size)])))
    
    plt.title(plot_titles[i])
    plt.legend(fontsize="10", loc="center left", bbox_to_anchor=(1, 0.5))
    
    if i==8 or i==1:
        plt.ylim((0, 10000))
        #plt.xlim((90, 200))

#plt.show()
print("Time taken to execute linear_sim once: ", run_time)
print("Estimated time to run 48,128 simulations: ", run_time*48128, "seconds, or ", run_time*48128/3600, "hours")