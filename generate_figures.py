"""
centralized figure generation file. Also acts as code documentation.
"""

# user modules
import sim_V2b as sv2
import project_library as PL
from M3_beta import beta_model_3

# system modules
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib
from pathlib import Path
import pandas as pd



#### set up LaTeX for plots.
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
#matplotlib.rcParams['font.size'] = 12

preamble = (r'\usepackage{amsmath}'
            r'\usepackage{siunitx}'
            r'\usepackage{bm}'
            r'\newcommand{\ve}{\varepsilon}')

matplotlib.rcParams['text.latex.preamble'] = preamble
fontsize = 12



def example():

    # Run(s)
    # ----------------------------------------------------------

    read_from_hyper = True
    read_from_param = True
    read_from_init = True

    save_init_states = False
    save_parameters = False
    save_hyperparams = False

    run_number = 'chronic_1'
    read_from_hyper_loc = Path.cwd() / 'presets' / 'hyper_preset_{}.csv'.format(run_number)
    read_from_param_loc = Path.cwd() / 'presets' / 'parameter_preset_{}.csv'.format(run_number)
    read_from_init_loc = Path.cwd() / 'presets' / 'init_val_preset_{}.csv'.format(run_number)
    

    hyp = sv2.get_hyp(read_from_hyper_loc,read_from_hyper,
                      run_number,save_hyperparams)
    init_dict = sv2.get_init_states(read_from_init_loc,read_from_init,
                                    run_number,save_init_states)
    parameters = sv2.get_parameters(read_from_param_loc,read_from_param,
                                    run_number,save_parameters)

    print(hyp)

    
    delta_t = hyp['delta_t']
    t_final = hyp['t_final']
    delayed_infection = hyp['delayed_infection']
    t_infection_start = hyp['t_infection_start']
    t_infection_end = hyp['t_infection_end']
    
    nosocomial_size = hyp['nosocomial_size']
    nosocomial_start = hyp['nosocomial_start']
    nosocomial_end = hyp['nosocomial_end']    

    HSPC_boost = hyp['HSPC_boost']
    HSPC_default = hyp['HSPC_default']
    HSPC_increment = hyp['HSPC_increment']
    HSPC_boost_time = hyp['HSPC_boost_time']

    SCSF_boost = hyp['SCSF_boost']
    SCSF_default = hyp['SCSF_default']
    SCSF_increment = hyp['SCSF_increment']
    SCSF_boost_time = hyp['SCSF_boost_time']


    
    init_state = [
        init_dict['HQ'],
        init_dict['HM'],
        init_dict['N'],
        init_dict['P'],
        init_dict['A'],
        init_dict['SCSF'],
        init_dict['K'],
        init_dict['Q'],
        init_dict['S'],
        init_dict['U'],
        init_dict['MDSC'],
        init_dict['MF']
    ]


    #data = pd.read_csv('your_file.csv')
    #print(data)

    path_increment = hyp['path_increment']
    path_size_default = hyp['path_size_default']

    outputs = []
    outputs_t = []

    for i in range(10):
        
        stim_times = [100]
        #stim_sizes = [path_size_default+i*path_increment]
        stim_sizes = [100+i*path_increment]

        t,data = PL.lin_sim_scipy(beta_model_3, parameters, init_state, t_final,
                                  delta_t, stim_times, stim_sizes)

        I = PL.calculate_I(data[3], data[4], data[6], data[2], parameters['theta_N'],
                           parameters['theta_K'], parameters['k'])

        data_I = np.concatenate([t.reshape(1,-1),data,I.reshape(1,-1)])

        outputs.append(data_I)


    # Plot
    # ----------------------------------------------------------

    fig,axs = plt.subplots(3,4,figsize=(8,7))
    axs1 = np.asarray(axs).flatten()

    titles = [
        '$H_Q(t)$',
        '$H_M(t)$',
        '$N(t)$',
        '$P(t)$',
        '$A(t)$',
        '$SCSF(t)$',
        '$K(t)$',
        '$Q(t)$',
        '$S(t)$',
        '$U(t)$',
        '$MDSC(t)$',
        '$MF(t)$', 
        '$I(t)$'
    ]



    for i in range(10):
        data_I = outputs[i]
        for j, out in enumerate(data_I[1:4]):
            axs1[j].plot(data_I[0], out, c=str(0.75*i/10),lw=1)
            axs1[j].set_title(titles[j])
    
    plt.tight_layout()
    
    
    return fig


def generate_figure(function, args, filenames, dpi=200):

    fig = function(*args)

    if type(filenames) is list:
        for name in filenames:
            fig.savefig(name,dpi=dpi)
    else:
        fig.savefig(filenames,dpi=dpi)


def main():

    # create figs directory if it doesn't exist
    if not(os.path.isdir('figs')):
        os.mkdir('figs')
    
    figures = [
        (example,[],['figs/example.pdf']),
        
    ]
    
    for fig in figures:
        generate_figure(*fig)


if __name__ == "__main__":
    main()
