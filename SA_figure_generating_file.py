import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# this is a template file, all that needs to be changed between output and order is the file names in step 1

# ========== 1: Load data from csv's into useable format ==========
time_array = np.arange(51, 75)

t_1 = pd.read_csv('total_Si_N_2_51.csv', delimiter='\t')
t_2 = pd.read_csv('total_Si_N_2_52.csv', delimiter='\t')
t_3 = pd.read_csv('total_Si_N_2_53.csv', delimiter='\t')
t_4 = pd.read_csv('total_Si_N_2_54.csv', delimiter='\t')
t_5 = pd.read_csv('total_Si_N_2_55.csv', delimiter='\t')
t_6 = pd.read_csv('total_Si_N_2_56.csv', delimiter='\t')
t_7 = pd.read_csv('total_Si_N_2_57.csv', delimiter='\t')
t_8 = pd.read_csv('total_Si_N_2_58.csv', delimiter='\t')
t_9 = pd.read_csv('total_Si_N_2_59.csv', delimiter='\t')
t_10 = pd.read_csv('total_Si_N_2_60.csv', delimiter='\t')
t_11 = pd.read_csv('total_Si_N_2_61.csv', delimiter='\t')
t_12 = pd.read_csv('total_Si_N_2_62.csv', delimiter='\t')
t_13 = pd.read_csv('total_Si_N_2_63.csv', delimiter='\t')
t_14 = pd.read_csv('total_Si_N_2_64.csv', delimiter='\t')
t_15 = pd.read_csv('total_Si_N_2_65.csv', delimiter='\t')
t_16 = pd.read_csv('total_Si_N_2_66.csv', delimiter='\t')
t_17 = pd.read_csv('total_Si_N_2_67.csv', delimiter='\t')
t_18 = pd.read_csv('total_Si_N_2_68.csv', delimiter='\t')
t_19 = pd.read_csv('total_Si_N_2_69.csv', delimiter='\t')
t_20 = pd.read_csv('total_Si_N_2_70.csv', delimiter='\t')
t_21 = pd.read_csv('total_Si_N_2_71.csv', delimiter='\t')
t_22 = pd.read_csv('total_Si_N_2_72.csv', delimiter='\t')
t_23 = pd.read_csv('total_Si_N_2_73.csv', delimiter='\t')
t_24 = pd.read_csv('total_Si_N_2_74.csv', delimiter='\t')

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

moderate_cutoff = np.zeros(24)+0.1
high_cutoff = np.zeros(24)+0.3

si_arrays = [gn_si_arr, knq_si_arr, kns_si_arr, ktn_si_arr, w_si_arr, p_crit_si_arr,
             spq_si_arr, sph_si_arr, sps_si_arr, sas_si_arr, sah_si_arr, sau_si_arr,
             Is_si_arr, As_si_arr, Ic_si_arr, Ac_si_arr, y_si_arr, ds_si_arr,
             dp_si_arr, da_si_arr, dq_si_arr, du_si_arr, H_init_si_arr, N_init_si_arr]

conf_arrays = [gn_conf_arr, knq_conf_arr, kns_conf_arr, ktn_conf_arr, w_conf_arr, p_crit_conf_arr,
             spq_conf_arr, sph_conf_arr, sps_conf_arr, sas_conf_arr, sah_conf_arr, sau_conf_arr,
             Is_conf_arr, As_conf_arr, Ic_conf_arr, Ac_conf_arr, y_conf_arr, ds_conf_arr,
             dp_conf_arr, da_conf_arr, dq_conf_arr, du_conf_arr, H_init_conf_arr, N_init_conf_arr]

titles = ['$g_{n}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$k_{nq}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$k_{ns}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$k_{tn}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$\omega$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$P_{crit}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$S_{pq}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$S_{ph}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$S_{ps}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$S_{as}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$S_{ah}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$S_{au}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$I_{S}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$A_{S}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$I_{c}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$A_{C}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$\gamma$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$d_{s}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$d_{p}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$d_{a}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$d_{q}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$d_{u}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$H_{Init}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)',
          '$N_{Init}$ Sensitivity Index During First 24 Hours of an Infection ($N(t)$)']

filenames = ['gn_SI_2', 'knq_SI_2', 'kns_SI_2', 'ktn_SI_2', 'w_SI_2', 'p_crit_SI_2',
             'spq_SI_2', 'sph_SI_2', 'sps_SI_2', 'sas_SI_2', 'sah_SI_2', 'sau_SI_2',
             'Is_SI_2', 'As_SI_2', 'Ic_SI_2', 'Ac_SI_2', 'y_SI_2', 'ds_SI_2',
             'dp_SI_2', 'da_SI_2', 'dq_SI_2', 'du_SI_2', 'H_init_SI_2', 'N_init_SI_2']

for i in range(24):

    fig, axs = plt.subplots()

    axs.plot(time_array, si_arrays[i], 'k')  # SI indices
    axs.plot(time_array, moderate_cutoff, 'y--', label='Mod. Influential lower bound')    # Lower bound on range to be considered moderately influential
    axs.plot(time_array, high_cutoff, 'g--', label='Highly Influential lower bound')        # Lower bound on range to be considered highly influential

    axs.fill_between(
    time_array, np.zeros(24)+si_arrays[i]-conf_arrays[i], np.zeros(24)+si_arrays[i]+conf_arrays[i], color='k', alpha=.15, label='$95\%$ Confidence Interval')

    axs.set_xlabel('Time ($t$)')
    axs.set_ylabel('Sensitivity Index')
    axs.set_ylim(ymin=0, ymax=1)
    axs.set_title(titles[i])
    axs.legend()

    fig.savefig(filenames[i])

    plt.close(fig)


print("Figures compiled succesfully!")
