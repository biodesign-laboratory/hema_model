import numpy as np
import matplotlib.pyplot as plt

from project_library import linear_sim_hybrid
from project_library import linear_sim_smooth

runs = 10
TFinal = 200
delta_t = 0.01
t = np.arange(0, TFinal, delta_t)

params_hybrid = [    
    
    1.2,    # g_N
    5,      # k_nq
    2,      # k_ns
    3,      # k_tn
    0.0005,  # w; UNIQUE to hybrid
    10000,   # p_crit; UNIQUE to hybrid
    0.5,      # S_PQ
    0.2,      # S_PH
    0.3,      # S_PS
    0.2,      # S_AS
    0.1,    # S_AH
    0.7,      # S_AU
    8000,   # immune_start; UNIQUE to hybrid
    6000,   # active_start; UNIQUE to hybrid
    20000,  # immune_crit; UNIQUE to hybrid
    15000,  # active_crit; UNIQUE to hybrid
    0.0005, # y; UNIQUE to hybrid
    0.3,    # d_s
    0.95,   # d_p
    0.8,    # d_a
    0.7,    # d_q
    0.85    # d_u

]

params_smooth = [

    1.2,    # g_n
    5,      # k_nq
    2,      # K_ns
    3,      # k_sn = k_tn
    5,      # S_PQ
    1,      # S_PH
    2,      # S_PS
    1,      # S_AS
    0.5,    # S_AH
    7,      # S_AU
    0.25,    # d_s
    0.95,   # d_p
    0.8,    # d_a
    0.3,    # d_q
    0.6,   # d_u
    1000,   # N_half; UNIQUE to smooth
    0.35,   # R_MAX; UNIQUE to smooth
    0.1,    # R_MIN; UNIQUE to smooth
    0.45,   # D_MAX; UNIQUE to smooth
    0,      # D_MIN; UNIQUE to smooth
    10000,  # R_CRIT; UNIQUE to smooth
    22000,  # D_CRIT; UNIQUE to smooth
    2000,   # S_half; UNIQUE to smooth
    2000,  # U_CRIT; UNIQUE to smooth
    500,    # Q_CRIT; UNIQUE to smooth
    -1500,   # Q_SHIFT
    -1500,   # U_SHIFT

]

init_values_hybrid = [

    1000,   # H_0
    0,      # N_0
    0,      # P_0
    0,      # A_0
    0,      # T_0 (combined total of all leukocytes)
    0,      # a_0 (percentage of total cells that are active)
    0.95,   # b_0 (percentage of total cells that are active)
    0.05,   # e_0 (percentage of total cells that are active)
    1       # E* initial
]

init_values_smooth = [

    1000,   # H_0
    0,      # N_0
    0,      # Q_0
    0,      # S_0
    0,      # U_0
    0,      # P_0
    0       # A_0
]

path_t = [100]  
path_size = np.zeros(5)

for i in range(len(path_size)):
    
    path_size[i] = 200 + 100*i
    
outputs_hybrid = np.zeros((5, 9, int(TFinal / delta_t)))
outputs_smooth = np.zeros((5, 9, int(TFinal / delta_t)))

for i in range(5):
    
    outputs_hybrid[i] = linear_sim_hybrid(init_values_hybrid, params_hybrid, delta_t, TFinal, path_t, [path_size[i]])
    outputs_smooth[i] = linear_sim_smooth(params_smooth, init_values_smooth, delta_t, TFinal, path_t, [path_size[i]])

graph_titles = [

    '$H(t)$',
    '$N(t)$',
    '$S(t)$',
    '$Q(t)$',
    '$U(t)$',
    '$P(t)$',
    '$A(t)$',
    '$T(t)$',
    '$I(t)$'

]

# fig, axs = plt.subplots(3, 3)

for i in range(9):
    
    plt.figure(i)
    for j in range(5):
        
        plt.plot(t, outputs_hybrid[j][i], label='hybrid')
        plt.plot(t, outputs_smooth[j][i], label='smooth')
        plt.title(graph_titles[i])
        plt.legend()
        plt.ylim(0, 3000)
        
        if(i == 3):
            plt.ylim(0, 500)


plt.show()
            
