import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import os
from scipy.integrate import solve_ivp


def reynolds_hema_ODE(t, y, p, stim_time=None):
    """
    ARGS:
    - p: dict containing relevant parameter values
    - y: array containing relevant state values [P(t), N(t), D(t), CA(t), M(t), I(t)]
    - t: to make compatible with solve_ivp
    --------------
    Returns arr [dP, dN, dD, dCA, dM, dI] of derivatives for each variable
    """
    k_pg = p["k_pg"]
    p_infty = p["p_infty"]
    k_pm = p["k_pm"]
    s_m = p['s_m']
    mu_m = p["mu_m"]
    k_mp = p["k_mp"]
    k_pn = p["k_pn"]
    mu_n = p["mu_n"]
    k_dn = p["k_dn"]
    x_dn = p["x_dn"]
    mu_d = p["mu_d"]
    s_c = p["s_c"]
    k_cn = p["k_cn"]
    k_cnd = p["k_cnd"]
    alpha_MS1 = p["alpha_MS1"]          # alpha parameter for M in S1
    alpha_M = p["alpha_M"]              # alpha parameter for M in C_A
    mu_c = p["mu_c"]
    alpha_D = p["alpha_D"]
    alpha_C = p["alpha_C"]
    mu_md = p["mu_md"]
    s_i = p["s_i"]
    mu_nr = p["mu_nr"]
    s_in = p["s_in"]
    C_infty = p["C_infty"]
    k_nn = p["k_nn"]
    k_np = p["k_np"]
    k_nd = p["k_nd"]
    
    P_t = y[0]
    N_t = y[1]
    D_t = y[2]
    CA_t = y[3]
    M_t = y[4]
    I_t = y[5]

    # relevant functions/quantities below
    def f(x):
        return x / (1 + (CA_t / C_infty)**2)
    
    def f_u(V, x):
        return V**6 / (x**6 + V**6)

    def f_d(V, x, n):
        return x**n / (x**n + V**n)
    
    R = f(k_nn*N_t+k_np*P_t+k_nd*D_t)

    S_1 = alpha_MS1*M_t + alpha_D*D_t + alpha_C*CA_t

    NR = (s_in*I_t)/(R + mu_nr)

    # derivates below

    dPdt = k_pg*P_t*(1 - (P_t/p_infty)) - (k_pm*s_m*P_t)/(mu_m + k_mp*P_t) - k_pn*f(N_t)*P_t    # Pathogens

    dNdt = R*NR - mu_n*N_t      # Neutrophils

    dDdt = k_dn*f_u(N_t, x_dn) - mu_d*D_t   # Damage

    dCAdt = s_c + (k_cn*f(N_t+k_cnd*D_t+alpha_M*M_t))/(1+f(N_t+k_cnd*D_t+alpha_M*M_t)) - mu_c*CA_t      # Anti-inflammatory mediators

    dMdt = S_1*I_t - mu_md*M_t       # MDSCs

    dIdt = s_i - mu_nr*I_t - s_in*I_t - S_1*I_t       # Immature Leukocytes

    return [dPdt, dNdt, dDdt, dCAdt, dMdt, dIdt]

def merge_figures_grid(nRow, nCol, img_width, img_height, exp_num, order, o_names, p_names, num=-1, directory=os.path.dirname(os.path.abspath(__file__))):
    '''
    ** template code, some lines will need to be modified as needed **
    Args:
    nRow: Number of rows in resulting image
    nCol: Number of columns in resulting image
    img_width: Image width of each individual image
    img_height: Image height of each individual image
    exp_num: Experiment number
    order: Sensitivity analysis specific argument, specifies which order graphs are being combined (first, second, or total)
    o_names: Array containing the desired output names
    p_names: Array containing the desired parameter names
    num: 
        -1 : Will not append any number to end of filename, only one merged figure created;
        else attached to end of filename; necessary if creating more than one merged figure for each parameter

    ====================
    Outputs:
    No output; Saves resulting image file, split into an nRow x nCol grid pattern composed of each individual image, into the
    same file location as the source code that called this function
    '''
    path_titles = p_names

    relative_path_arr = np.empty(len(o_names), dtype=f'<U256')    # will contain paths to the images in the same order that the output names appear in o_names arg

    script_dir = directory

    for param in path_titles:

        for i, out_name in enumerate(o_names):
            
            filepath = os.path.join(script_dir, f'Experiment_{exp_num}', f'{out_name}_out', f'{order}_figs', f'{param}_SI_{exp_num}_{out_name}.png')
            # print(f'{param} for {out_name} exists: ' + str(os.path.exists(filepath)))
            relative_path_arr[i] = os.path.join(script_dir, f'Experiment_{exp_num}', f'{out_name}_out', f'{order}_figs', f'{param}_SI_{exp_num}_{out_name}.png')

        # Create a new blank image with a white background
        collage_width = nCol * img_width
        collage_height = nRow * img_height
        collage_image = Image.new('RGB', (collage_width, collage_height), 'white')

        # Paste each image into the collage
        for i, image_file in enumerate(relative_path_arr):
            # print(f"Opening image: {image_file}")
            if os.path.exists(image_file):
                try:
                    img = Image.open(image_file)
                    # print(f"Image size before resizing: {img.size}")
                    img = img.resize((img_width, img_height), Image.Resampling.LANCZOS)
                    # print(f"Image size after resizing: {img.size}")
                    x = (i % nCol) * img_width
                    y = (i // nCol) * img_height
                    # print(f"Pasting image at coordinates: ({x}, {y})")
                    collage_image.paste(img, (x, y))
                except Exception as e:
                    print(f"Error opening or pasting image: {e}")
            else:
                print(f"Image file does not exist: {image_file}")

    # Save the collage image
        filepath = os.path.join(script_dir, f'Experiment_{exp_num}', f'{order}_order_merge')
        if not os.path.exists(filepath):
            os.makedirs(filepath)

        if num != -1:
            collage_image.save(os.path.join(filepath, f'{param}_merge_{order}.png'))        # weird behavior going on here, still worked somehow, must check
        else:
            collage_image.save(os.path.join(filepath, f'{param}_merge_{order}.png'))

def csv_to_figure(o_names, o_name_laTex, p_names, p_names_laTex, nTimesteps, init_time, nDatapoints, order, exp_num, filepath='default'):

    '''
    o_names: array containing strings of output names used for file paths
    o_names_laTex: same as above but used for formatting matplotlib plot titles in laTex
    p_names: array containing strings of input names used for file paths
    p_names_laTex: same as above but used for formatting matplotlib plot titles in laTex
    init_time: initial time
    nDatapoints: number of data types (i.e. columns) to include from .csv
    order: case-specific to sensitivity analysis
    exp_num: trial number
    filepath: OPTIONAL, string containing parent directory of all relevant data folders, 'default' arg sets path to relative path

    This function is specifically meant to be used with sobol script
    ==========================
    Does not return anything, generates time-series figures from given .csv files and stores them in predefined folder locations
    '''
    script_dir = ''
    if filepath=='default':

        script_dir = os.path.dirname(os.path.abspath(__file__))     # file path is relative to source code path

    else:
        script_dir = filepath


    nParam = len(p_names)

    # ========== 1: Load data from csv's into useable format ==========

    master_df = np.zeros((len(o_names), nTimesteps, nParam, nDatapoints))
    # master_df: output -> SIs for all params sorted by time -> SIs sorted by param e.g. to access the sensitivity index of parameter 'z' in the timestep 'y' for output 'x', the index would be df[x, y, z, 1]

    for i, str in enumerate(o_names):

        SIs_per_timestep = np.zeros((nTimesteps, nParam, nDatapoints))

        for t in np.arange(init_time, init_time + nTimesteps):

            output_timestep_Si_data = pd.read_csv(os.path.join(script_dir, f'{str}_out', f'{order}_order', f'{order}_Si_{str}_{exp_num}_{t}.csv'), delimiter='\t')  # dataframe with 24 rows, 3 columns where rows=parameter, column0=param_name, column1=Si_index, and column2=SI_conf
            output_timestep_Si_data = output_timestep_Si_data.iloc[:, [1, 2]]
            SIs_per_timestep[t-init_time] = output_timestep_Si_data.to_numpy()
            
        master_df[i] = SIs_per_timestep

    # =============== 2. Plot data =============================

    moderate_cutoff = np.zeros(nTimesteps)+0.1
    high_cutoff = np.zeros(nTimesteps)+0.3

    titles = np.empty(len(p_names), dtype=f'<U256')
    for i, p_name in enumerate(p_names_laTex):

        titles[i] = f'{p_name}'

    
    filenames = np.empty(len(p_names), dtype=f'<U256')

    for p, p_name in enumerate(p_names):
        
        filenames[p] = f'{p_name}_SI_{exp_num}_'

    for i, out_name in enumerate(o_names):

        for p, param_name in enumerate(p_names):

            fig, axs = plt.subplots()

            axs.plot(np.arange(init_time, init_time+nTimesteps), master_df[i, :, p, 0], 'k')  # SI indices
            axs.plot(np.arange(init_time, init_time+nTimesteps), moderate_cutoff, 'y--', label='Mod. Influential lower bound')    # Lower bound on range to be considered moderately influential
            axs.plot(np.arange(init_time, init_time+nTimesteps), high_cutoff, 'g--', label='Highly Influential lower bound')        # Lower bound on range to be considered highly influential
            axs.plot(np.arange(init_time, init_time+nTimesteps), np.zeros(nTimesteps), 'k--')

            axs.fill_between(np.arange(init_time, init_time+nTimesteps), master_df[i, :, p, 0]-master_df[i, :, p, 1], master_df[i, :, p, 0]+master_df[i, :, p, 1], color='k', alpha=.15, label='$95\%$ Confidence Interval')
            # axs.axvline(x=50, color='r', linestyle='--', label='Pathogen Input Timestep')    # x must be changed depending on pathogen input time, add arg to toggle this on or off

            axs.set_xlabel('Time ($t$)')
            axs.set_ylabel('Sensitivity Index')
            axs.set_ylim(ymin=-0.05, ymax=1.05)
            axs.set_title(f'${titles[p]}$' + " Sensitivity Index for Output " + f'${o_name_laTex[i]}(t)$')
            axs.legend(loc = 'upper right')

            path = os.path.join(script_dir, f'{out_name}_out', f'{order}_figs')

            if not os.path.exists(path):
                os.makedirs(path)

            fig.savefig(os.path.join(path, f'{param_name}_SI_{exp_num}_{out_name}.png'))

            plt.close(fig)

def event_function(t, y, parameters, stim_time):
    """
    trigger event when t-delta_time = 0
    space is just a dummy variable. the number of parameters
    after t,y need to match the terms in the "args" keyword arg
    in the solve_ivp call.
    """
    
    return t - stim_time

def lin_sim_scipy(ODE_eq, parameters, y0, tf, dt, stim_times=[], stim_sizes=[], force_timesteps=False):
    
    '''
    lin_sim for scipy
    assume single delta function stimuli for now.
    Next: heaviside inputs and multiple stimuli with index choice.
    
    ARGS:
    ODE_eq : Function argument, first-order derivatives of system output
    parameters : dict, must be compatible with ODE_eq
    init_y : Initial state of the system, init_y=(y_1(0), y_2(0), ..., y_n(0))
    t_final : Final timestep to reach
    delta_t : Timestep size

    stim_time : time(s) of delta-function
    stim_size : stimulus size(s)
    
    force_timesteps : forces solve_ivp to return output along all timepoints using t_eval; used for sensitivity analysis

    OUTPUT:
    solution
    '''


    event_function.terminal = True
    t = np.arange(0,tf,dt)
    data = np.zeros([1,len(y0)])
    times = np.zeros(1)

    for i,stim_time in enumerate(stim_times):
        
        if force_timesteps:
            out = solve_ivp(ODE_eq,(t[0],t[-1]),y0,
                        args=(parameters,stim_time),
                        events=event_function, t_eval=t)
        else:
            out = solve_ivp(ODE_eq,(t[0],t[-1]),y0,
                        args=(parameters,stim_time),
                        events=event_function)

        # save values just before first pert.
        data = np.concatenate([data,out.y.T])
        times = np.concatenate([times,out.t])

        # update inputs for next iteration
        t = np.arange(out.t[-1]+dt,tf,dt)
        y0 = out.y[:,-1]
        y0[2] += stim_sizes[i]

    if force_timesteps:
        out = solve_ivp(ODE_eq,(t[0],t[-1]),y0,
                    args=(parameters,stim_time),
                    events=event_function, t_eval=t)
    else:
        out = solve_ivp(ODE_eq,(t[0],t[-1]),y0,
                        args=(parameters,stim_time),
                        events=event_function)
    data = np.concatenate([data,out.y.T])
    times = np.concatenate([times,out.t])

    return times[1:],data.T[:,1:]

    

    """
    # start terminal sim
    event_function.terminal = True

    t = np.arange(0,tf,dt)    
    
    out = solve_ivp(ODE_eq,(t[0],t[-1]),y0,
                    args=(parameters,stim_times[0]),
                    events=event_function)

    # terminate on first stim time

    # restart sim

    y2 = out.y[:,-1]
    y2[2] += stim_sizes[0]
    t2 = np.arange(out.t[-1]+dt,tf,dt)
    print(out.t[-1],tf)

    out2 = solve_ivp(ODE_eq,(t2[0],t2[-1]),y2,
                     args=(parameters,stim_times[0]),
                     method='LSODA')

    # terminate on next stim time

    data = np.concatenate([out.y.T,out2.y.T])
    t = np.concatenate([out.t,out2.t])
    print(out.t[-1],out2.t[0],out2.t[-1])

    return t,data.T
            
    """
