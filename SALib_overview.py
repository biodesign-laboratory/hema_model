import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from SALib.test_functions import Ishigami
from project_library import linear_sim
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

param_values_IHD = sobol_sample.sample(problem_IHD, 1024)
print("Shape of the generated sample: ", param_values_IHD.shape)        # sanity check
print("First few samples:")
print(param_values_IHD[:3])

init_values = np.hstack((param_values_IHD[:, -2:], np.tile(np.array([4000, 0, 0.9, 0.1, 1800, 300, 1]), (param_values_IHD.shape[0], 1))))
print("First few initial value sets: ")         # sanity check
print(init_values[3,:])
print("Shape of init_value array: ", init_values.shape)

IHD_out = np.zeros((param_values_IHD.shape[0], 9))

print("! Computations starting now !")
start_time = time.time()
for i, X in enumerate(param_values_IHD):

    IHD_out[i] = linear_sim(init_values[i], X[:23], 0.01, 200, [100], [5000])[:, -1]        # for each function in the output data, grab only the last time-series data point

end_time = time.time()
elapsed_time = end_time - start_time
print("Outputs generated for all input samples. Elapsed time: ", elapsed_time)

# save input samples and output to external files to avoid having to calculate outputs for first experiment
np.savetxt('input_samples.txt', param_values_IHD)
np.savetxt('output_from_samples.txt', IHD_out)
np.savetxt('elapsed_time.txt', np.array([elapsed_time]))

# these will be saved in the SA folder ^^

Si = sobol_analyze.analyze(problem_IHD, IHD_out)

total_Si, first_Si, second_Si = Si.to_df()
total_Si.to_csv('total_Si.to_csv', sep='\t')
first_Si.to_csv('first_Si.to_csv', sep='\t')
second_Si.to_csv('second_Si.to_csv', sep='\t')
