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
timestep_size = 0.1
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
        N_vals[i] = N_0 + i*100 + 3500
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


# ---------------- plots: -------------------
num_graphs = 9
# -- Base plots for each output from 0 <= t <= TFinal --
for i in range(num_graphs):
    
    plt.figure(i)
    for output in outputs:
        
        plt.plot(t, output[i], label="N_0: " + str(output[1][0]))
    
    plt.title(plot_titles[i])
    plt.legend(fontsize="10", loc="center left", bbox_to_anchor=(1, 0.5))

# -- Plots focused on pathogen insults


'''z = 10*(runs-1)+8   # +8 because we have 8 function outputs

#fig,axs = plt.subplots(8)
plt_figs = np.zeros(runs)       # left off here

for i in range(z):      # avoids double for loop

    if i % 10 >= 8:
        continue
    
    else:

        fig = plt.figure(i//10)
        plt.plot(t, outputs[i//10][i%10])
        plt.title(plot_titles[i%10])
        plt.legend(fontsize="8", loc="center left", bbox_to_anchor=(1, 0.5))
        plt_figs[i//10] = fig'''




''' ------------------------------- construction above, safe code below --------------
plt.title(plot_titles[i])
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
plt.ylim(0, 2000)

plt.figure(1)
plt.plot(t, output_0[0])
plt.plot(t, output_1[0])
plt.plot(t, output_2[0])
plt.plot(t, output_3[0])
plt.plot(t, output_4[0])
plt.plot(t, output_5[0])
plt.plot(t, output_6[0])
plt.plot(t, output_7[0])
plt.plot(t, output_8[0])
plt.plot(t, output_9[0])
plt.plot(t, output_10[0])
#plt.plot(t, [1000]*200, "k--")
plt.title('HSPC\'s')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
plt.ylim(0, 2000)


plt.figure(2)
plt.plot(t, output_0[1], label="N_0: " + str(INIT_VALUES_0[1]))
plt.plot(t, output_1[1], label="N_0: " + str(INIT_VALUES_1[1]))
plt.plot(t, output_2[1], label="N_0: " + str(INIT_VALUES_2[1]))
plt.plot(t, output_3[1], label="N_0: " + str(INIT_VALUES_3[1]))
plt.plot(t, output_4[1], label="N_0: " + str(INIT_VALUES_4[1]))
plt.plot(t, output_5[1], label="N_0: " + str(INIT_VALUES_5[1]))
plt.plot(t, output_6[1], label="N_0: " + str(INIT_VALUES_6[1]))
plt.plot(t, output_7[1], label="N_0: " + str(INIT_VALUES_7[1]))
plt.plot(t, output_8[1], label="N_0: " + str(INIT_VALUES_8[1]))
plt.plot(t, output_9[1], label="N_0: " + str(INIT_VALUES_9[1]))
plt.plot(t, output_10[1],label= "N_0: " + str(INIT_VALUES_10[1]))
#plt.plot(t, [N_0]*200, "k--")
plt.title('Pathogens (95<=t<=150)')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
plt.xlim([95, 150])
plt.ylim([0, 10000])
plt.legend(fontsize="8")

plt.figure(3)
plt.plot(t, output_0[2])
plt.plot(t, output_1[2])
plt.plot(t, output_2[2])
plt.plot(t, output_3[2])
plt.plot(t, output_4[2])
plt.plot(t, output_5[2])
plt.plot(t, output_6[2])
plt.plot(t, output_7[2])
plt.plot(t, output_8[2])
plt.plot(t, output_9[2])
plt.plot(t, output_10[2])
#plt.plot(t, np.ones(200), "k--")
plt.title('Pro-inflammatory cytokines')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(4)
plt.plot(t, output_0[3])
plt.plot(t, output_1[3])
plt.plot(t, output_2[3])
plt.plot(t, output_3[3])
plt.plot(t, output_4[3])
plt.plot(t, output_5[3])
plt.plot(t, output_6[3])
plt.plot(t, output_7[3])
plt.plot(t, output_8[3])
plt.plot(t, output_9[3])
plt.plot(t, output_10[3])
plt.title('Anti-inflammatory cytokines')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(4)
plt.plot(t, output_0[5], "r-")
plt.plot(t, output_0[7], "b-")
plt.plot(t, output_0[4], "g-")
plt.plot(t, output_0[6], "k-")
plt.title("Active / Total / Stable / Immuno Leukocytes")
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
plt.legend(["Active", "Total", "Stable", "Immuno-"], loc ="lower right")

plt.figure(5)
plt.plot(t, output_0[4])
plt.plot(t, output_1[4])
plt.plot(t, output_2[4])
plt.plot(t, output_3[4])
plt.plot(t, output_4[4])
plt.plot(t, output_5[4])
plt.plot(t, output_6[4])
plt.plot(t, output_7[4])
plt.plot(t, output_8[4])
plt.plot(t, output_9[4])
plt.plot(t, output_10[4])
plt.title('Stable Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(6)
plt.plot(t, output_0[5])
plt.plot(t, output_1[5])
plt.plot(t, output_2[5])
plt.plot(t, output_3[5])
plt.plot(t, output_4[5])
plt.plot(t, output_5[5])
plt.plot(t, output_6[5])
plt.plot(t, output_7[5])
plt.plot(t, output_8[5])
plt.plot(t, output_9[5])
plt.plot(t, output_10[5])
plt.title('Active Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(7)
plt.plot(t, output_0[6])
plt.plot(t, output_1[6])
plt.plot(t, output_2[6])
plt.plot(t, output_3[6])
plt.plot(t, output_4[6])
plt.plot(t, output_5[6])
plt.plot(t, output_6[6])
plt.plot(t, output_7[6])
plt.plot(t, output_8[6])
plt.plot(t, output_9[6])
plt.plot(t, output_10[6])
plt.title('Immuno-suppressive Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(8)
plt.plot(t, output_0[7])
plt.plot(t, output_1[7])
plt.plot(t, output_2[7])
plt.plot(t, output_3[7])
plt.plot(t, output_4[7])
plt.plot(t, output_5[7])
plt.plot(t, output_6[7])
plt.plot(t, output_7[7])
plt.plot(t, output_8[7])
plt.plot(t, output_9[7])
plt.plot(t, output_10[7])
plt.title('Total Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(9)
plt.plot(t, output_0[0], label="N_0: " + str(INIT_VALUES_0[1]))
plt.plot(t, output_1[0], label="N_0: " + str(INIT_VALUES_1[1]))
plt.plot(t, output_2[0], label="N_0: " + str(INIT_VALUES_2[1]))
plt.plot(t, output_3[0], label="N_0: " + str(INIT_VALUES_3[1]))
plt.plot(t, output_4[0], label="N_0: " + str(INIT_VALUES_4[1]))
plt.plot(t, output_5[0], label="N_0: " + str(INIT_VALUES_5[1]))
plt.plot(t, output_6[0], label="N_0: " + str(INIT_VALUES_6[1]))
plt.plot(t, output_7[0], label="N_0: " + str(INIT_VALUES_7[1]))
plt.plot(t, output_8[0], label="N_0: " + str(INIT_VALUES_8[1]))
plt.plot(t, output_9[0], label="N_0: " + str(INIT_VALUES_9[1]))
plt.plot(t, output_10[0],label= "N_0: " + str(INIT_VALUES_10[1]))
#plt.plot(t, [N_0]*200, "k--")
plt.title('HSPCs (0<=t<=10)')
plt.xlabel("t")
plt.ylabel("Density")
plt.xlim([0, 10])
plt.ylim([0, 1100])
plt.legend(fontsize="8")

plt.figure(10)
plt.plot(t, output_0[0], label="N_0: " + str(INIT_VALUES_0[1]))
plt.plot(t, output_1[0], label="N_0: " + str(INIT_VALUES_1[1]))
plt.plot(t, output_2[0], label="N_0: " + str(INIT_VALUES_2[1]))
plt.plot(t, output_3[0], label="N_0: " + str(INIT_VALUES_3[1]))
plt.plot(t, output_4[0], label="N_0: " + str(INIT_VALUES_4[1]))
plt.plot(t, output_5[0], label="N_0: " + str(INIT_VALUES_5[1]))
plt.plot(t, output_6[0], label="N_0: " + str(INIT_VALUES_6[1]))
plt.plot(t, output_7[0], label="N_0: " + str(INIT_VALUES_7[1]))
plt.plot(t, output_8[0], label="N_0: " + str(INIT_VALUES_8[1]))
plt.plot(t, output_9[0], label="N_0: " + str(INIT_VALUES_9[1]))
plt.plot(t, output_10[0],label= "N_0: " + str(INIT_VALUES_10[1]))
plt.title('HSPCs (95<=t<=150)')
plt.xlabel("t")
plt.ylabel("Density")
plt.xlim([95, 150])
plt.ylim([0, 1100])
plt.legend(fontsize="8", loc="center left", bbox_to_anchor=(1, 0.5))

plt.figure(11)
plt.plot(t, output_0[1], label="N_0: " + str(INIT_VALUES_0[1]))
plt.plot(t, output_1[1], label="N_0: " + str(INIT_VALUES_1[1]))
plt.plot(t, output_2[1], label="N_0: " + str(INIT_VALUES_2[1]))
plt.plot(t, output_3[1], label="N_0: " + str(INIT_VALUES_3[1]))
plt.plot(t, output_4[1], label="N_0: " + str(INIT_VALUES_4[1]))
plt.plot(t, output_5[1], label="N_0: " + str(INIT_VALUES_5[1]))
plt.plot(t, output_6[1], label="N_0: " + str(INIT_VALUES_6[1]))
plt.plot(t, output_7[1], label="N_0: " + str(INIT_VALUES_7[1]))
plt.plot(t, output_8[1], label="N_0: " + str(INIT_VALUES_8[1]))
plt.plot(t, output_9[1], label="N_0: " + str(INIT_VALUES_9[1]))
plt.plot(t, output_10[1],label= "N_0: " + str(INIT_VALUES_10[1]))
plt.title('Pathogens (0<=t<=15)')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
plt.xlim([0, 15])
plt.ylim([0, 6000])
plt.legend(fontsize="8")

plt.figure(12)
plt.plot(t, output_n_steady[1])
plt.title("Pathogens")
#plt.ylim([n_steady-1000, n_steady+1000])

print(n_steady)

'''