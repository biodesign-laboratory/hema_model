# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
This is the main function that calculates the values for each individual function from the derivatives using linear approximation
"""
# This function.py is placed in the SA folder; it has been updated (as of 6/26) with the following changes:
#   - Removed 8 obsolete parameters
#   - Fixed N_inf, S_A, and S_N 


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import os


def linear_sim(init_values, rates, timestep_size, TFinal, path_repeat_t, path_repeat_size):
    '''

    Parameters
    ----------
    init_values : Initial values for each function; f(T=0)
    
    rates : parameter values 
    
    timestep_size : Size of each timestep; as t_count -> 0, runtime and accuracy increase
    
    TFinal : last t value (an integer)
    
    path_repeat_t : time(s) to repeat a pathogen infection; must be array
    
    path_repeat_size : strength of pathogen force at each t; must be same size as path_repeat_t

    Returns
    -------
    output = ndarray of shape (#, t_count); output[0] = ndarray of timesteps, 
    output[i] denotes the py.list containing the values of function f_i for each t in output[0]
    
    returns 0 if error

    '''
    #--------- 0. Verifying params work ---------------
    
    rates_len = len(rates) >= 22
    path_repeat_lens = len(path_repeat_t) == len(path_repeat_size)
    bools = [rates_len, path_repeat_lens]
    
    for i in bools:
    
        if i == False:
            return 0
    
    #--------- 1. Initiliazing all values and parameters
    
    timesteps = np.arange(0, TFinal, timestep_size)
    
    H_0 = init_values[0]        # initial HSPC value
    N_0 = init_values[1]        # initial pathogen value
    T_0 = init_values[2]        # initial combined value of all leukocytes
    a_0 = init_values[3]        # initial alpha value
    b_0 = init_values[4]        # initial beta value
    e_0 = init_values[5]        # initial epsilon value
    P_0 = init_values[6]        # initial pro-inflammatory cytokines value
    A_0 = init_values[7]        # initial anti-inflammatory cytokines value
    E_star_0 = init_values[8]       # initial E_star value (added 1/12)
    
    H = H_0
    N = N_0
    T = T_0
    a = a_0
    b = b_0
    e = e_0
    P = P_0
    A = A_0
    E_star = E_star_0
    
    Q = a * T
    S = b * T
    U = e * T
    
    H_output = [H]
    N_output = [N]
    T_output = [T]
    Q_output = [Q]
    S_output = [S]
    U_output = [U]
    P_output = [P]
    A_output = [A]
    T_output = [T]
    
    g_N = rates[0]                  # pathogen growth coefficient per hour
    k_nq = rates[1]                 # kill rate of pthogens by active WBCs
    k_ns = rates[2]                 # kill rate of pthogens by stable WBCs
    k_tn = rates[3]                 # kill rate of stable WBCs by stable leukocytes
    w = rates[4]                    # how quickly does E(t) -> 0 (as a function of inflammation I(t)); considering fixing
    P_crit = rates[5]               # value of inflammation at which 
    S_a = 1                         # fixed (SA)
    S_n = 1                         # fixed (SA)
    N_inf = 2*10**7                 # fixed (SA)
    S_PQ = rates[6]                 # secretion rate of pro by active WBCs
    S_PH = rates[7]                 # secretion rate of pro by HSPCs
    S_PS = rates[8]                 # secretion rate of pro by stable WBCs
    S_AS = rates[9]                 # secretion rate of anti by stable WBCs
    S_AH = rates[10]                # secretion rate of anti by HSPCs
    S_AU = rates[11]                # secretion rate of anti by immuno-suppressive WBCs
    Immune_start = rates[12]        # value of inflammation at which stable->immuno-suppresive WBC conversion begins
    Active_start = rates[13]        # value of inflammation at which stable->active WBC conversion begins
    Immune_crit = rates[14]         # value of inflammation at which stable->immuno-suppresive WBC reaches its max value
    Active_crit = rates[15]         # value of inflammation at which stable->active WBC conversion reaches its max value
    y = rates[16]                   # modulates the pathogen killing rate of stable WBCs by the density of stable WBCs
    d_s = rates[17]                 # natural decay rate of stable WBCs
    d_p = rates[18]                 # natural decay rate of pro-inflammatory cytokines
    d_a = rates[19]                 # natural decay rate of anti-inflammatory cytokines
    d_q = rates[20]                 # natural decay rate of active WBCs
    d_u = rates[21]                 # natural decay rate of immuno_suppressive WBCs
    
    count = 1
    
    I_output = [P_0+S_n*N_0-S_a*A_0]
    #----------- 2. Calculating derivative values for each t > t_0
    
    for x in timesteps[1:]:
        
            
        #---------- mechanism for repeating pathogenic insult ---------------------
        
        for i in range(len(path_repeat_t)):
            if count == path_repeat_t[i]/timestep_size:
                N += path_repeat_size[i]
        
        count += 1
        
        #---------- 2a. Calculating individual terms ---------------
        "Function for E (exhaustion - dH/dt)"
        
        if (P + S_n * N - S_a * A) < P_crit:
            E = 1
        elif (P + S_n * N - S_a * a) >= P_crit:
            try:
                E = 2 - (2/(1+np.exp(-1 * y * (P + S_n * N - S_a * A - P_crit))))
            except FloatingPointError:
                E = 0
            
        if E < 0:
            E = 0
            
        "Function for R (self renewal - dH/dt)"
        
        if P + S_n * N - S_a * A <= 0:
            R = 0
        elif P + S_n * N - S_a * A <= P_crit:
            R = (0.35 * ((P + S_n * N - S_a * A)/P_crit)) * H
        elif P + S_n * N - S_a * A > P_crit:
            R = 0.35 * H
            
        "Function for D (differentiation - dH/dt & dS/dt)"
        
        if (P + S_n * N - S_a * A <= 0):
            D = 0
        elif ((P + S_n * N - S_a * A) <= P_crit):
            D = (0.35*((P + S_n * N - S_a * A)/P_crit))*H
        elif ((P + S_n * N - S_a * a) > P_crit):
            D = 0.35*H
            
        "Determining E_star for H_stable"
        
        if P + S_n*N - S_a*A <= 1.5*P_crit or E_star <= 0.1:
            dE_star = 0
            
        else:
            dE_star = -0.0005*(P + S_n*N - S_a*A - 1.5*P_crit)
            
        "Linear approximation for E_star(t)"
        
        E_star = E_star + timestep_size*dE_star
        if E_star < 0.1:
            E_star = 0.1
        
        
        "BM Niche Renewal contribution (dH/dt)"
        
        if E_star * H_0 > H:
            H_stable = 0.05 * H * (1-(H/(E_star * H_0)))
        
        else:
            H_stable = 0
            
        "Function for A_l (dS/dt - stable leukocyte kill rate by pathogens)"
        
        try:
            A_l = (k_tn * N) * ((-2/(1 + np.exp(w * S))) + 1)
        except FloatingPointError:
            A_l = k_tn * N

        "Function for mu_SP"
        
        if P + S_n*N - S_a*A <= Active_start:
            mu_SP = 0
            
        elif P + S_n*N - S_a*A <= Active_crit:
            mu_SP = 0.5 * ((P + S_n*N - S_a*A - Active_start) / Active_crit) * S
            
        else:
            mu_SP = 0.5 * S

            
        "Function for mu_SA"
        
        if P + S_n*N - S_a*A <= Immune_start:
            mu_SA = 0
            
        elif P + S_n*N - S_a*A <= Immune_crit:
            mu_SA = 0.2 * ((P + S_n*N - S_a*A - Immune_start) / Immune_crit) * S
            
        else:
            mu_SA = 0.2 * S
        
        
        "Function for D_P (dP/dt - pro-cytokine decay term"

        D_P = d_p * P
        
        "Function for D_A (dA/dt - anti-cytokine decay term"
        
        if A > 0:
            D_A = d_a * A
        else:
            D_A = 0
        
        
        "Function for D_S (S decay)"
        
        D_S = d_s * S
        
        "Function for D_Q (Q decay)"
        
        D_Q = d_q * Q
        
        "Function for D_U (U decay)"
        
        D_U = d_u * U
        
        
        #---------- 2b. Calculating derivatives ---------------------
        
        dHdt = E*(R + H_stable) - D                                     # HSPC derivative
        
        dNdt = (g_N*N - ((k_nq * Q) + (k_ns * S)) * ((-2/(1 + np.exp(w * S))) + 1)) * (1-(N/N_inf))        # Pathogen derivative
        
        dPdt = (S_PS * S) + (S_PQ * Q) + (S_PH * H) - D_P     # Pro-inflammatory derivative
        
        dAdt = (S_AU * U) + (S_AS * S) + (S_AH * H) - D_A     # Anti-inflammatory derivative
        
        dSdt = D - A_l - D_S - mu_SA - mu_SP                            # Stable Leukocytes derivative
        
        dQdt = mu_SP - D_Q                                  # Active Leukocytes derivative
        
        dUdt = mu_SA - D_U                                  # Immuno-suppressive derivative
        
        #--------- 2c. Diagnostics --------------------------------
        
        #fill with code for checking variable values as needed
        
        #--------- 3. Updating lists and functions with linear approximation -------------
        
        H = timestep_size*dHdt + H        # control flow may affect function values
        if H < 0:
            H = 0
        
        N = timestep_size*dNdt + N
        if N < 0:
            N = 0
        
        S = timestep_size*dSdt + S
        if S < 0:
            S = 0
        
        Q = timestep_size*dQdt + Q
        if Q < 0:
            Q = 0
        
        U = timestep_size*dUdt + U
        if U < 0:
            U = 0
                
        P = timestep_size*dPdt + P
        if P < 0:
            P = 0
        
        A = timestep_size*dAdt + A
        if A < 0:
            A = 0
        
        T = Q + U + S
        
        I = P+S_n*N-S_a*A
        
        H_output.append(H)
        N_output.append(N)
        P_output.append(P)
        A_output.append(A)
        S_output.append(S)
        Q_output.append(Q)
        U_output.append(U)
        T_output.append(T)
        I_output.append(I)
        
        # the condition to set =0 if X<0 is due to the discrete nature of linear approximation; if X=0 is a fixed point then the
        # smooth functions should, in theory, never go negative
    #------------- 4. Output --------------
    
    output = np.array([H_output, N_output, P_output, A_output, S_output, Q_output, U_output, T_output, I_output])
    
    return output
    
        
def merge_figures_grid(nRow, nCol, img_width, img_height, file_str_arr):
    '''
    ** template code, some lines will need to be modified as needed **
    Args:
    nRow: Number of rows in resulting image
    nCol: Number of columns in resulting image
    img_width: Image width of each individual image
    img_height: Image height of each individual image
    file_str_arr: Array containing relevant file title information for each file to merge
    ====================
    Outputs:
    No output; Saves resulting image file, split into an nRow x nCol grid pattern composed of each individual image, into the
    same file location as the source code that called this function
    '''
    path_titles = ['Ac', 'As', 'da', 'dp', 'dq', 'ds', 'du', 'gn', 
               'H_init', 'Ic', 'Is', 'knq', 'kns', 'ktn', 'N_init', 'p_crit',
               'sah', 'sas', 'sau', 'sph', 'spq', 'sps', 'w', 'y']

    for param in path_titles:

        script_dir = os.path.dirname(os.path.abspath(__file__))
        relative_path_H = os.path.join(script_dir, 'H_out', 'first_order', 'figs', f'{param}_SI_2.png')
        relative_path_I = os.path.join(script_dir, 'I_out', 'first_order', 'figs', f'{param}_SI_2.png')
        relative_path_N = os.path.join(script_dir, 'N_out', 'first_order', 'figs', f'{param}_SI_2.png')
        relative_path_Q = os.path.join(script_dir, 'Q_out', 'first_order', 'figs', f'{param}_SI_2.png')
        relative_path_S = os.path.join(script_dir, 'S_out', 'first_order', 'figs', f'{param}_SI_2.png')
        relative_path_U = os.path.join(script_dir, 'U_out', 'first_order', 'figs', f'{param}_SI_2.png')


        # Define the grid size and image size
        grid_columns = nCol
        grid_rows = nRow
        image_width = img_width  # Width of each individual image
        image_height = img_height  # Height of each individual image

        # Create a new blank image with a white background
        collage_width = grid_columns * image_width
        collage_height = grid_rows * image_height
        collage_image = Image.new('RGB', (collage_width, collage_height), 'white')

        # List of image file paths
        image_files = [
            relative_path_H, 
            relative_path_I, 
            relative_path_N,
            relative_path_Q, 
            relative_path_S, 
            relative_path_U
        ]

        # Paste each image into the collage
        for i, image_file in enumerate(image_files):
            img = Image.open(image_file)
            img = img.resize((image_width, image_height), Image.Resampling.LANCZOS)
            x = (i % grid_columns) * image_width
            y = (i // grid_columns) * image_height
            collage_image.paste(img, (x, y))

        # Save the collage image
        collage_image.save(f'{param}_merge_first.png')


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# this is a template file, the code's purpose is to generate time-series line graphs and save them as .png's from given sensitivity analysis output .csv data
# for file paths, this file is expected to be in the same directory as the folders containing output data

def csv_to_figure(output_names, param_names, nTimesteps, init_time, nDatapoints, order, exp_num, filepath='no_path'):
    '''
    filepath: OPTIONAL, string containing parent directory of all relevant data folders, 'no_path' arg sets path to relative path
    output_names: array containing strings of output names
    param_names: array containing strings of input names
    init_time: initial time
    nDatapoints: number of data types (i.e. columns) to include from .csv
    order: case-specific to sensitivity analysis
    exp_num: trial number

    WARNING: order='second' not fully implemented
    ==========================
    No output, generates time-series figures from given .csv files
    '''
    script_dir = ''
    if filepath=='no_path':

        script_dir = os.path.dirname(os.path.abspath(__file__))     # file path is relative to source code path

    else:
        script_dir = filepath


    nParam = len(param_names)

    # ========== 1: Load data from csv's into useable format ==========

    master_df = np.zeros((len(output_names), nTimesteps, nParam, nDatapoints))
    # master_df: output -> SIs for all params sorted by time -> SIs sorted by param e.g. to access the sensitivity index of paramater 'z' in the timestep 'y' for output 'x', the index would be df[x, y, z, 1]

    for i, str in enumerate(output_names):

        SIs_per_timestep = np.zeros((nTimesteps, nParam, nDatapoints))

        for t in np.arange(init_time, init_time + nTimesteps):

            output_timestep_Si_data = pd.read_csv(os.path.join(script_dir, f'{str}_out', f'{order}_order', f'{order}_Si_{str}_{exp_num}_{t}.csv'), delimiter='\t')  # dataframe with 24 rows, 3 columns where rows=parameter, column0=param_name, column1=Si_index, and column2=SI_conf
            output_timestep_Si_data = output_timestep_Si_data.iloc[:, [1, 2]]
            SIs_per_timestep[t-51] = output_timestep_Si_data.to_numpy()
            
        master_df[i] = SIs_per_timestep

    # =============== 2. Plot data =============================

    moderate_cutoff = np.zeros(24)+0.1
    high_cutoff = np.zeros(24)+0.3

    # outputs_in_laTex = ['($H(t)$)', '($I(t)$)', '($N(t)$)', '($Q(t)$)', '($S(t)$)', '($U(t)$)']

    outputs_in_laTex = np.empty(len(output_names), dtype='str')
    for i, out_name in enumerate(output_names):

        outputs_in_laTex[i] = f'${out_name}(t)'

    '''titles = np.empty(len(param_names), dtype='str'):
    for i, p_name in enumerate(param_names):

        titles[i] = f'$'''

    titles = ['$g_{n}$ Sensitivity Index During First 24 Hours of an Infection',
            '$k_{nq}$ Sensitivity Index During First 24 Hours of an Infection',
            '$k_{ns}$ Sensitivity Index During First 24 Hours of an Infection',
            '$k_{tn}$ Sensitivity Index During First 24 Hours of an Infection',
            '$\omega$ Sensitivity Index During First 24 Hours of an Infection',
            '$P_{crit}$ Sensitivity Index During First 24 Hours of an Infection',
            '$S_{pq}$ Sensitivity Index During First 24 Hours of an Infection',
            '$S_{ph}$ Sensitivity Index During First 24 Hours of an Infection',
            '$S_{ps}$ Sensitivity Index During First 24 Hours of an Infection',
            '$S_{as}$ Sensitivity Index During First 24 Hours of an Infection',
            '$S_{ah}$ Sensitivity Index During First 24 Hours of an Infection',
            '$S_{au}$ Sensitivity Index During First 24 Hours of an Infection',
            '$I_{S}$ Sensitivity Index During First 24 Hours of an Infection',
            '$A_{S}$ Sensitivity Index During First 24 Hours of an Infection',
            '$I_{C}$ Sensitivity Index During First 24 Hours of an Infection',
            '$A_{C}$ Sensitivity Index During First 24 Hours of an Infection',
            '$\gamma$ Sensitivity Index During First 24 Hours of an Infection',
            '$d_{s}$ Sensitivity Index During First 24 Hours of an Infection',
            '$d_{p}$ Sensitivity Index During First 24 Hours of an Infection',
            '$d_{a}$ Sensitivity Index During First 24 Hours of an Infection',
            '$d_{q}$ Sensitivity Index During First 24 Hours of an Infection',
            '$d_{u}$ Sensitivity Index During First 24 Hours of an Infection',
            '$H_{Init}$ Sensitivity Index During First 24 Hours of an Infection',
            '$N_{Init}$ Sensitivity Index During First 24 Hours of an Infection']

    filenames = [f'gn_SI_{exp_num}_', f'knq_SI_{exp_num}_', f'kns_SI_{exp_num}_', f'ktn_SI_{exp_num}_', f'w_SI_{exp_num}_', f'p_crit_SI_{exp_num}_',
                f'spq_SI_{exp_num}_', f'sph_SI_{exp_num}_', f'sps_SI_{exp_num}_', f'sas_SI_{exp_num}_', f'sah_SI_{exp_num}_', f'sau_SI_{exp_num}_',
                f'Is_SI_{exp_num}_', f'As_SI_{exp_num}_', f'Ic_SI_{exp_num}_', f'Ac_SI_{exp_num}_', f'y_SI_{exp_num}_', f'ds_SI_{exp_num}_',
                f'dp_SI_{exp_num}_', f'da_SI_{exp_num}_', f'dq_SI_{exp_num}_', f'du_SI_{exp_num}_', f'H_init_SI_{exp_num}_', f'N_init_SI_{exp_num}_']

    for i, out_name in enumerate(output_names):

        for p, param_name in enumerate(param_names):    # enumerate is unnecessary here, edit later

            fig, axs = plt.subplots()

            axs.plot(np.arange(init_time, init_time+nTimesteps), master_df[i, :, p, 0], 'k')  # SI indices
            axs.plot(np.arange(init_time, init_time+nTimesteps), moderate_cutoff, 'y--', label='Mod. Influential lower bound')    # Lower bound on range to be considered moderately influential
            axs.plot(np.arange(init_time, init_time+nTimesteps), high_cutoff, 'g--', label='Highly Influential lower bound')        # Lower bound on range to be considered highly influential

            axs.fill_between(
            np.arange(init_time, init_time+nTimesteps), master_df[i, :, p, 0]-master_df[i, :, p, 1], master_df[i, :, p, 0]+master_df[i, :, p, 1], color='k', alpha=.15, label='$95\%$ Confidence Interval')

            axs.set_xlabel('Time ($t$)')
            axs.set_ylabel('Sensitivity Index')
            axs.set_ylim(ymin=0, ymax=1)
            axs.set_title(titles[p] + " " + outputs_in_laTex[i])
            axs.legend()

            path = os.path.join(script_dir, f'{out_name}_out', f'{order}_figs')

            if not os.path.exists(path):
                os.makedirs(path)

            fig.savefig(os.path.join(path, f'{param_name}_SI_{exp_num}_{out_name}.png'))

            plt.close(fig)


    print("Figures compiled succesfully!")
