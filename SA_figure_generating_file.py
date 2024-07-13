import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# this is a template file, all that needs to be changed between output and order is the file names in step 1

# ========== 1: Load data from csv's into useable format ==========
time_array = np.arange(51, 75)

t_1 = pd.read_csv('total_Si_H_2_51.csv', delimiter='\t')
t_2 = pd.read_csv('total_Si_H_2_52.csv', delimiter='\t')
t_3 = pd.read_csv('total_Si_H_2_53.csv', delimiter='\t')
t_4 = pd.read_csv('total_Si_H_2_54.csv', delimiter='\t')
t_5 = pd.read_csv('total_Si_H_2_55.csv', delimiter='\t')
t_6 = pd.read_csv('total_Si_H_2_56.csv', delimiter='\t')
t_7 = pd.read_csv('total_Si_H_2_57.csv', delimiter='\t')
t_8 = pd.read_csv('total_Si_H_2_58.csv', delimiter='\t')
t_9 = pd.read_csv('total_Si_H_2_59.csv', delimiter='\t')
t_10 = pd.read_csv('total_Si_H_2_60.csv', delimiter='\t')
t_11 = pd.read_csv('total_Si_H_2_61.csv', delimiter='\t')
t_12 = pd.read_csv('total_Si_H_2_62.csv', delimiter='\t')
t_13 = pd.read_csv('total_Si_H_2_63.csv', delimiter='\t')
t_14 = pd.read_csv('total_Si_H_2_64.csv', delimiter='\t')
t_15 = pd.read_csv('total_Si_H_2_65.csv', delimiter='\t')
t_16 = pd.read_csv('total_Si_H_2_66.csv', delimiter='\t')
t_17 = pd.read_csv('total_Si_H_2_67.csv', delimiter='\t')
t_18 = pd.read_csv('total_Si_H_2_68.csv', delimiter='\t')
t_19 = pd.read_csv('total_Si_H_2_69.csv', delimiter='\t')
t_20 = pd.read_csv('total_Si_H_2_70.csv', delimiter='\t')
t_21 = pd.read_csv('total_Si_H_2_71.csv', delimiter='\t')
t_22 = pd.read_csv('total_Si_H_2_72.csv', delimiter='\t')
t_23 = pd.read_csv('total_Si_H_2_73.csv', delimiter='\t')
t_24 = pd.read_csv('total_Si_H_2_74.csv', delimiter='\t')

t_x_arr = [t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, t_10,
           t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, t_19,
           t_20, t_21, t_22, t_23, t_24,]   # array for easy indexing purposes

# =========== 2. Create the index and conf_int arrays for each parameter ==============

gn_si_arr = np.zeros(24)        # si_arr's contain the sensitivity indices for each parameter at each timestep
knq_si_arr = np.zeros(24)
kns_si_arr = np.zeros(24)
ktn_si_arr = np.zeros(24)
w_si_arr = np.zeros(24)
p_crit_si_arr = np.zeros(24)
spq_si_arr = np.zeros(24)
sph_si_arr = np.zeros(24)
sps_si_arr = np.zeros(24)
sas_si_arr = np.zeros(24)
sah_si_arr = np.zeros(24)
sau_si_arr = np.zeros(24)
Is_si_arr = np.zeros(24)
As_si_arr = np.zeros(24)
Ic_si_arr = np.zeros(24)
Ac_si_arr = np.zeros(24)
y_si_arr = np.zeros(24)
ds_si_arr = np.zeros(24)
dp_si_arr = np.zeros(24)
da_si_arr = np.zeros(24)
dq_si_arr = np.zeros(24)
du_si_arr = np.zeros(24)
H_init_si_arr = np.zeros(24)
N_init_si_arr = np.zeros(24)

gn_conf_arr = np.zeros(24)        # conf_arr's contain the confidence interval data for each parameter at each timestep
knq_conf_arr = np.zeros(24)
kns_conf_arr = np.zeros(24)
ktn_conf_arr = np.zeros(24)
w_conf_arr = np.zeros(24)
p_crit_conf_arr = np.zeros(24)
spq_conf_arr = np.zeros(24)
sph_conf_arr = np.zeros(24)
sps_conf_arr = np.zeros(24)
sas_conf_arr = np.zeros(24)
sah_conf_arr = np.zeros(24)
sau_conf_arr = np.zeros(24)
Is_conf_arr = np.zeros(24)
As_conf_arr = np.zeros(24)
Ic_conf_arr = np.zeros(24)
Ac_conf_arr = np.zeros(24)
y_conf_arr = np.zeros(24)
ds_conf_arr = np.zeros(24)
dp_conf_arr = np.zeros(24)
da_conf_arr = np.zeros(24)
dq_conf_arr = np.zeros(24)
du_conf_arr = np.zeros(24)
H_init_conf_arr = np.zeros(24)
N_init_conf_arr = np.zeros(24)

# ========== 3. Load data into the arrays ======================

for t in range(24):

    gn_si_arr[t] = t_x_arr[t].iloc[0, 1]
    knq_si_arr[t] = t_x_arr[t].iloc[1, 1]
    kns_si_arr[t] = t_x_arr[t].iloc[2, 1]
    ktn_si_arr[t] = t_x_arr[t].iloc[3, 1]
    w_si_arr[t] = t_x_arr[t].iloc[4, 1]
    p_crit_si_arr[t] = t_x_arr[t].iloc[5, 1]
    spq_si_arr[t] = t_x_arr[t].iloc[6, 1]
    sph_si_arr[t] = t_x_arr[t].iloc[7, 1]
    sps_si_arr[t] = t_x_arr[t].iloc[8, 1]
    sas_si_arr[t] = t_x_arr[t].iloc[9, 1]
    sah_si_arr[t] = t_x_arr[t].iloc[10, 1]
    sau_si_arr[t] = t_x_arr[t].iloc[11, 1]
    Is_si_arr[t] = t_x_arr[t].iloc[12, 1]
    As_si_arr[t] = t_x_arr[t].iloc[13, 1]
    Ic_si_arr[t] = t_x_arr[t].iloc[14, 1]
    Ac_si_arr[t] = t_x_arr[t].iloc[15, 1]
    y_si_arr[t] = t_x_arr[t].iloc[16, 1]
    ds_si_arr[t] = t_x_arr[t].iloc[17, 1]
    dp_si_arr[t] = t_x_arr[t].iloc[18, 1]
    da_si_arr[t] = t_x_arr[t].iloc[19, 1]
    dq_si_arr[t] = t_x_arr[t].iloc[20, 1]
    du_si_arr[t] = t_x_arr[t].iloc[21, 1]
    H_init_si_arr[t] = t_x_arr[t].iloc[22, 1]
    N_init_si_arr[t] = t_x_arr[t].iloc[23, 1]

    gn_conf_arr[t] = t_x_arr[t].iloc[0, 2]
    knq_conf_arr[t] = t_x_arr[t].iloc[1, 2]
    kns_conf_arr[t] = t_x_arr[t].iloc[2, 2]
    ktn_conf_arr[t] = t_x_arr[t].iloc[3, 2]
    w_conf_arr[t] = t_x_arr[t].iloc[4, 2]
    p_crit_conf_arr[t] = t_x_arr[t].iloc[5, 2]
    spq_conf_arr[t] = t_x_arr[t].iloc[6, 2]
    sph_conf_arr[t] = t_x_arr[t].iloc[7, 2]
    sps_conf_arr[t] = t_x_arr[t].iloc[8, 2]
    sas_conf_arr[t] = t_x_arr[t].iloc[9, 2]
    sah_conf_arr[t] = t_x_arr[t].iloc[10, 2]
    sau_conf_arr[t] = t_x_arr[t].iloc[11, 2]
    Is_conf_arr[t] = t_x_arr[t].iloc[12, 2]
    As_conf_arr[t] = t_x_arr[t].iloc[13, 2]
    Ic_conf_arr[t] = t_x_arr[t].iloc[14, 2]
    Ac_conf_arr[t] = t_x_arr[t].iloc[15, 2]
    y_conf_arr[t] = t_x_arr[t].iloc[16, 2]
    ds_conf_arr[t] = t_x_arr[t].iloc[17, 2]
    dp_conf_arr[t] = t_x_arr[t].iloc[18, 2]
    da_conf_arr[t] = t_x_arr[t].iloc[19, 2]
    dq_conf_arr[t] = t_x_arr[t].iloc[20, 2]
    du_conf_arr[t] = t_x_arr[t].iloc[21, 2]
    H_init_conf_arr[t] = t_x_arr[t].iloc[22, 2]
    N_init_conf_arr[t] = t_x_arr[t].iloc[23, 2]

# =============== 4. Plot data =============================

fig, axs = plt.subplots()

# plot g_n
axs.plot(time_array, gn_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('g_n_sensitivity_over_time.png')
plt.clf()

# plot k_nq
axs.plot(time_array, knq_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('k_nq_sensitivity_over_time.png')
plt.clf()

# plot k_ns
axs.plot(time_array, kns_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('k_ns_sensitivity_over_time.png')
plt.clf()

# plot k_tn
axs.plot(time_array, ktn_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('k_tn_sensitivity_over_time.png')
plt.clf()

# plot w
axs.plot(time_array, w_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('w_sensitivity_over_time.png')
plt.clf()

# plot p_crit
axs.plot(time_array, p_crit_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('p_crit_sensitivity_over_time.png')
plt.clf()

# plot s_pq
axs.plot(time_array, spq_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('s_pq_sensitivity_over_time.png')
plt.clf()

# plot sph
axs.plot(time_array, sph_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('s_ph_sensitivity_over_time.png')
plt.clf()

# plot sps
axs.plot(time_array, sps_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('s_ps_sensitivity_over_time.png')
plt.clf()

# plot sas
axs.plot(time_array, sas_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('s_as_sensitivity_over_time.png')
plt.clf()

# plot s_ah
axs.plot(time_array, sah_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('s_ah_sensitivity_over_time.png')
plt.clf()

# plot s_au
axs.plot(time_array, sau_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('s_au_sensitivity_over_time.png')
plt.clf()

# plot I_s
axs.plot(time_array, Is_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('I_s_sensitivity_over_time.png')
plt.clf()

# plot A_s
axs.plot(time_array, As_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('A_s_sensitivity_over_time.png')
plt.clf()

# plot I_c
axs.plot(time_array, Ic_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('I_c_sensitivity_over_time.png')
plt.clf()

# plot A_c
axs.plot(time_array, Ac_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('A_c_sensitivity_over_time.png')
plt.clf()

# plot y
axs.plot(time_array, y_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('y_sensitivity_over_time.png')
plt.clf()

# plot d_s
axs.plot(time_array, ds_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('d_s_sensitivity_over_time.png')
plt.clf()

# plot d_p
axs.plot(time_array, dp_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('d_p_sensitivity_over_time.png')
plt.clf()

# plot d_a
axs.plot(time_array, da_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('d_a_sensitivity_over_time.png')
plt.clf()

# plot d_q
axs.plot(time_array, dq_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('d_q_sensitivity_over_time.png')
plt.clf()

# plot d_u
axs.plot(time_array, du_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('d_u_sensitivity_over_time.png')
plt.clf()

# plot H_init
axs.plot(time_array, H_init_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('H_init_sensitivity_over_time.png')
plt.clf()

# plot N_init
axs.plot(time_array, N_init_si_arr)
axs.set_ylim( ymin=0, ymax=1)
plt.savefig('N_init_sensitivity_over_time.png')
plt.clf()

print("Figures compiled succesfully!")
