import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from SALib.test_functions import Ishigami
from project_library import linear_sim

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
    'num_vars': 22,
    'names': ['g_n', 'k_nq', 'k_ns', 'k_tn', 'w', 'p_crit', 's_pq', 's_ph', 's_ps', 's_as', 's_ah', 's_au', 
               'I_s', 'A_s', 'I_c', 'A_c', 'y', 'd_s', 'd_p', 'd_a', 'd_q', 'd_u'],
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
               [0.5, 1.5]]
}

init_values = [1000, 0, 4000, 0, 0.9, 0.1, 1800, 300, 1] # chosen mostly arbitrarily

param_values_IHD = sobol_sample.sample(problem_IHD, 1024)
print("Shape of the generated sample: ", param_values_IHD.shape)
print("First few samples:")
print(param_values_IHD[:3])

IHD_out = np.zeros((param_values_IHD.shape[0], 9))   # for this test, I'm only evaluating sensitivity on H(t)

for i, X in enumerate(param_values_IHD):
    IHD_out[i] = linear_sim(init_values, X, 0.01, 200, [100], [5000])[:,-1]         # only grab time-series output of H(t)


