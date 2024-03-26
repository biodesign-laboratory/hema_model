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
k_nq = 2
k_ns = 0.5
k_tn = 0.33
w = 0.0005
P_crit = 7000
A_crit = 1000
t_half_leukocytes = 7
t_half_P = 4.1375
t_half_A = 7
t_double_P = 15
gamma = 0.0000001
S_a = 1
S_n = 1
N_inf = 20000000
S_PQ = 2
S_PH = 0.33
S_PS = 0.02
S_AS = 0.005
S_AH = 0.005
S_AU = 2
theta_ps = 10_000_000
theta_ar = 10_000_000
theta_AS = 10_000_000
theta_UP = 10_000_000
Immune_start = 2000
Active_start = 3000
Immune_crit = 7000
Active_crit = 5000
y = 0.000005
d_s = 0.1
d_p = 0.75
d_a = 0.75
d_q = 0.8
d_u = 0.9
z=1.55

RATES = [g_N, K_PS, K_AQ, K_AS, K_PU, k_nq, k_ns, k_tn, w, P_crit, A_crit, S_a, S_n, N_inf, S_PQ, S_PH, S_PS, S_AS, S_AH, S_AU, 
         theta_ps, theta_ar, theta_AS, theta_UP, Immune_start, Active_start, Immune_crit, Active_crit, y, d_s, 
         d_p, d_a, d_q, d_u, z]

TFinal = 200
timestep_size = 0.05
t = np.arange(0, TFinal, timestep_size)


H_0 = 5000
N_0 = 0
P_0 = 1800
A_0 = 300
T_0 = 4000
a_0 = 0.00
b_0 = 0.973
e_0 = 0.027
E_star_0 = 1


N_vals = np.zeros(runs)

plot_label_nums = np.zeros(runs)

path_repeat_ts = [50, 100, 150]         # this should not be changed
path_repeat_arr = []

for i in range(runs):
        
    if i == 0:      # no pathogen insult, steady-state control run
        path_repeat_size = [0, 0, 0]
        path_repeat_arr.append([path_repeat_ts, path_repeat_size])
    elif 1 <= i <= 3:   # low pathogen insult, return to homeostasis
        path_repeat_size = [0, 6500+500*i, 0]
        path_repeat_arr.append([path_repeat_ts, path_repeat_size])
    else:
        path_repeat_size = [0, 10500+150*i, 0]
        path_repeat_arr.append([path_repeat_ts, path_repeat_size])

init_vals = np.eye(runs, 9)

for i in range(runs):
    init_vals[i] = np.array([H_0, N_0, T_0, a_0, b_0, e_0, P_0, A_0, E_star_0])



outputs = np.zeros((runs, 10, int(TFinal//timestep_size+1)))


for i in range(runs):
    
    outputs[i]=simulate_TEST(init_vals[i], RATES, timestep_size, TFinal, path_repeat_arr[i][0], path_repeat_arr[i][1])



plot_titles = [ 'HSPC\'s' , 'External Force', 'Pro-inflammatory cytokines', 'Anti-inflammatory cytokines',
                 'Stable WBCs', 'Activated WBCs', 'Immunosuppressive WBCs', 'Total WBCs', 'Total Inflammation', 'dE_star']


# ---------------- plots: -------------------
num_graphs = 9
# -- Base plots for each output from 0 <= t <= TFinal --
for i in range(num_graphs):
    
    plt.figure(i)
    for output in outputs:
        
        plt.plot(t, output[i], label="$N_{100}$: " + str(int(output[1][(int(100/timestep_size))])))
    
    plt.title(plot_titles[i])
    plt.legend(fontsize="10", loc="center left", bbox_to_anchor=(1, 0.5))
    
    if i==8 or i==1:
        
        A = P_crit+np.zeros(int(TFinal/timestep_size))
        B = 1.5*A
        plt.plot(t, A, 'k--', label="$P_{crit}$")
        plt.plot(t, B, 'r--', label="$1.5*P_{crit}$")
        plt.ylim((0, 15000))
        plt.xlim((90, 200))
        
    if i == 9:
        #plt.ylim((1300, 1600))
        plt.xlim((90, 150))

    plt.legend(fontsize="10", loc="center left", bbox_to_anchor=(1, 0.5))
    
# for diagnosing issues with simulation, output[i] for all i >=10 represents a
# debugging array. num_graphs=9 will give only the model outputs needed
