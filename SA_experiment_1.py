import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from SALib.test_functions import Ishigami

from project_library import linear_sim_smooth as linear_sim
import time
import pandas

# code block directly below is following SALib tutorial
'''problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-1*np.pi, np.pi],
               [-1*np.pi, np.pi],
               [-1*np.pi, np.pi]]
}       # dictionary defining the inputs to perform SA on

param_values = saltelli.sample(problem, 1024)   # generates 1024*(2*3+2) input samples

Y = Ishigami.evaluate(param_values)     # evaluates the Ishigami function on generated input samples

Si = sobol.analyze(problem, Y)'''
# Si is a py dictionary with keys "S1", "S2", "ST", "S1_conf", "S2_conf", "ST_conf"
# the "_conf" keys store the corresponding confidence intervals (95%)

'''print("First-order:", Si['S1'])
print("Second-order:", Si['S2'])
print("Total-order:", Si['ST'])'''


# -----------------------------------------------------------------------------------------------------------------
# code block directly below is using linear_sim function to evaluate sensitivity on the IHD model 
problem_IHD = {
    'num_vars': 24,
    'names': ['g_n', 'k_nq', 'k_ns', 'k_tn', 'w', 'p_crit', 's_pq', 's_ph', 's_ps', 's_as', 's_ah', 's_au', 
               'I_s', 'A_s', 'I_c', 'A_c', 'y', 'd_s', 'd_p', 'd_a', 'd_q', 'd_u', 'H_init', 'N_init'],
    'bounds': [[0.5, 2.5],
               [3, 10],
               [1, 3],
               [0.5, 2],
               [0.00005, 0.005],
               [3000, 10000],
               [2, 5],
               [0.25, 1],
               [0.5, 2],
               [0.5, 2],
               [0.25, 1],
               [2, 5],
               [1000, 5000],
               [2000, 6000],
               [4000, 10000],
               [3000, 8000],
               [0.00005, 0.005],
               [0.05, 0.25],
               [0.75, 1.5],
               [0.5,1.25],
               [0.75, 2],
               [0.5, 1.5],
               [500, 5000],
               [0, 20000]]
}

# init_values = [1000, 0, 4000, 0, 0.9, 0.1, 1800, 300, 1] # chosen mostly arbitrarily

param_values_IHD = sobol_sample.sample(problem_IHD, 2000)
print("Shape of the generated sample: ", param_values_IHD.shape)
print("First few samples:")
print(param_values_IHD[:3])

init_values = np.hstack((param_values_IHD[:, -2:-1], np.tile(np.array([0, 4000, 0, 0.9, 0.1, 1800, 300, 1]), (param_values_IHD.shape[0], 1))))
print("First few initial value sets: ")         # sanity check
print(init_values[3,:])
print("Shape of init_value array: ", init_values.shape)

IHD_out = np.zeros((param_values_IHD.shape[0], 9))

print("! Computations starting now !")
start_time = time.time()
for i, X in enumerate(param_values_IHD):

    IHD_out[i] = linear_sim(init_values[i], X[:23], 0.01, 100, [50], [X[-1]])[:, -1]        # for each function in the output data, grab only the last time-series data point

end_time = time.time()
elapsed_time = end_time - start_time
print("Outputs generated for all input samples. Elapsed time: ", elapsed_time)


SI_H = sobol_analyze.analyze(problem_IHD, IHD_out[:,0].flatten())     # sensitivity indices for H(t) (HSPCs)
SI_N = sobol_analyze.analyze(problem_IHD, IHD_out[:,1].flatten())     # sensitivity indices for N(t) (Pathogens)
SI_S = sobol_analyze.analyze(problem_IHD, IHD_out[:,4].flatten())     # sensitivity indices for S(t) (Stable cells)
SI_Q = sobol_analyze.analyze(problem_IHD, IHD_out[:,5].flatten())     # sensitivity indices for Q(t) (Active cells)
SI_U = sobol_analyze.analyze(problem_IHD, IHD_out[:,6].flatten())     # sensitivity indices for U(t) (Immuno-suppressive cells)
SI_I = sobol_analyze.analyze(problem_IHD, IHD_out[:,8].flatten())     # sensitivity indices for I(t) (Net Inflammation)


total_Si_H, first_Si_H, second_Si_H = SI_H.to_df()
total_Si_N, first_Si_N, second_Si_N = SI_N.to_df()
total_Si_S, first_Si_S, second_Si_S = SI_S.to_df()
total_Si_Q, first_Si_Q, second_Si_Q = SI_Q.to_df()
total_Si_U, first_Si_U, second_Si_U = SI_U.to_df()
total_Si_I, first_Si_I, second_Si_I = SI_I.to_df()

total_Si_H.to_csv('total_Si_H.csv', sep='\t')
first_Si_H.to_csv('first_Si_H.csv', sep='\t')
second_Si_H.to_csv('second_Si_H.csv', sep='\t')

total_Si_N.to_csv('total_Si_N.csv', sep='\t')
first_Si_N.to_csv('first_Si_N.csv', sep='\t')
second_Si_N.to_csv('second_Si_N.csv', sep='\t')

total_Si_S.to_csv('total_Si_S.csv', sep='\t')
first_Si_S.to_csv('first_Si_S.csv', sep='\t')
second_Si_S.to_csv('second_Si_S.csv', sep='\t')

total_Si_Q.to_csv('total_Si_Q.csv', sep='\t')
first_Si_Q.to_csv('first_Si_Q.csv', sep='\t')
second_Si_Q.to_csv('second_Si_Q.csv', sep='\t')

total_Si_U.to_csv('total_Si_U.csv', sep='\t')
first_Si_U.to_csv('first_Si_U.csv', sep='\t')
second_Si_U.to_csv('second_Si_U.csv', sep='\t')

total_Si_I.to_csv('total_Si_I.csv', sep='\t')
first_Si_I.to_csv('first_Si_I.csv', sep='\t')
second_Si_I.to_csv('second_Si_I.csv', sep='\t')
