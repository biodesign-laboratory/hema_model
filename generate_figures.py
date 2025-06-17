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
import matplotlib as mpl
import time

import string # for labels
labels = list(string.ascii_uppercase)



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

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Run(s)
    # ----------------------------------------------------------

    plot_aseptic = True        # False = plot infection

    read_from_hyper = True
    read_from_param = True
    read_from_init = True

    save_init_states = False
    save_parameters = False
    save_hyperparams = False

    '''if not Path.exists(Path.cwd() / 'figure_presets'):
        Path.mkdir(Path.cwd() / 'figure_presets')'''

    run_number = 'aseptic_tau_Q'
    read_from_hyper_loc = Path.cwd() / 'presets' / 'hyper_preset_{}.csv'.format(run_number)
    read_from_param_loc = Path.cwd() / 'presets' / 'parameter_preset_{}.csv'.format(run_number)
    read_from_init_loc = Path.cwd() / 'presets' / 'init_val_preset_{}.csv'.format(run_number)
    
    hyp = sv2.get_hyp(read_from_hyper_loc,read_from_hyper,
                      run_number,save_hyperparams, plot_aseptic)
    
    init_dict = sv2.get_init_states(read_from_init_loc,read_from_init,
                                    run_number,save_init_states)
    parameters = sv2.get_parameters(read_from_param_loc,read_from_param,
                                    run_number,save_parameters)


    runs = hyp['runs']
    delta_t = hyp['delta_t']
    t_final = hyp['t_final']

    if not plot_aseptic:
        t_infection_start = hyp['t_infection_start']
        nosocomial_size = hyp['nosocomial_size']
        nosocomial_start = hyp['nosocomial_start']
        path_increment = hyp['path_increment']
        path_size_default = hyp['path_size_default']

        if nosocomial_start <= t_infection_start:
            raise Exception("Nosocomial must occur after the initial infection")
    else:
        parameter_to_increment = hyp['parameter_to_increment']
        increment_by = hyp['increment_by']

    
    init_state = [
        init_dict['HST'],
        init_dict['MPP'],
        init_dict['N'],
        init_dict['P'],
        init_dict['A'],
        init_dict['SCSF'],
        init_dict['K'],
        init_dict['Q'],
        init_dict['S'],
        init_dict['U'],
        init_dict['MDSC'],
        init_dict['EN']
    ]

    if plot_aseptic:

        param_arr = [parameters]

        for i in range(1, runs):
            x = parameters.copy()
            x[parameter_to_increment] = parameters[parameter_to_increment] + increment_by*i
            param_arr.append(x)

    else:

        pathogen_sizes = np.arange(path_size_default,
                                path_size_default + runs*path_increment,
                                path_increment)

    outputs = []
        

    for i in range(runs):

        if not plot_aseptic:        # for infection graphs
            stim_times = [t_infection_start, nosocomial_start]
            stim_sizes = [pathogen_sizes[i], nosocomial_size]

            start = time.time()

            t,data = PL.lin_sim_scipy(beta_model_3, parameters, init_state, t_final,
                                    delta_t, stim_times, stim_sizes)
            
            end = time.time()
            print(f"Run {i+1} successfully calculated, time taken: {end-start}s")

            I = PL.calculate_I(data[3], data[4], data[6], data[2], parameters['theta_N'],
                            parameters['theta_K'], parameters['k'])

            data_I = np.concatenate([t.reshape(1,-1),data,I.reshape(1,-1)])

            outputs.append(data_I)
        
        else:       # for aseptic graphs

            stim_times=[100]
            stim_size=[0]

            start = time.time()

            t,data = PL.lin_sim_scipy(beta_model_3, param_arr[i], init_state, t_final,
                                    delta_t, stim_times, stim_size)
            
            end = time.time()
            print(f"Run {i+1} successfully calculated, time taken: {end-start}s")

            I = PL.calculate_I(data[3], data[4], data[6], data[2], param_arr[i]['theta_N'],
                            param_arr[i]['theta_K'], param_arr[i]['k'])

            data_I = np.concatenate([t.reshape(1,-1),data,I.reshape(1,-1)])

            outputs.append(data_I)



    # Plot
    # ----------------------------------------------------------

    fig,axs = plt.subplots(3,3,figsize=(8,6))
    axs1 = np.asarray(axs).flatten()

    output_indices_to_plot = [1, 2, 3, 4, 5, 7, 8, 9, 10]       # choose which outputs to plot

    titles = [
        '$H_{ST}(t)$',
        '$MPP(t)$',
        '$N(t)$',
        '$P(t)$',
        '$A(t)$',
        '$SCSF(t)$',
        '$K(t)$',
        '$Q(t)$',
        '$S(t)$',
        '$U(t)$',
        '$MDSC(t)$',
        '$EN(t)$', 
        '$I(t)$'
    ]


    cmap = mpl.colormaps['viridis']
    colors = cmap(np.linspace(0, .9, runs))

    for i in range(runs):
        data_I = outputs[i]
        
        for j, index in enumerate(output_indices_to_plot):
            
            out = data_I[index]
            if j == 0 and not plot_aseptic:
                label = str(pathogen_sizes[i])

            else:
                label = ''       # exact pathogen loads may not be important for figures, commented out for now
                
            axs1[j].plot(data_I[0], out, c=colors[i],lw=1,label=label)

            if i == 0:
                axs1[j].set_title(labels[j],loc='left')
                axs1[j].set_xlabel('$t$')
                axs1[j].set_ylabel(titles[index-1])

                if not plot_aseptic:
                    for k in range(2):
                        
                        if k == 1:
                            if nosocomial_size != 0:
                                axs1[j].axvline(x=stim_times[k],color='orange',lw=1,ls='--')
                        else:
                            axs1[j].axvline(x=stim_times[k],color='orange',lw=1,ls='--')

    # make some plots on log scale
    axs1[2].set_yscale('log')
    axs1[2].set_ylim((1.0,1e8))

    #axs1[0].legend(fontsize='small',labelspacing=0.1)
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
        (example,[],[f'figs/debug.pdf']),
        (example,[],[f'figs/debug.png'])
    ]
    
    for fig in figures:
        generate_figure(*fig)


if __name__ == "__main__":
    main()
