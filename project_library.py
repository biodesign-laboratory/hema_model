import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import os


def linear_sim_hybrid(init_values, rates, timestep_size, TFinal, path_repeat_t, path_repeat_size):
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
    S_a = rates[22]                 # strength of anti-inflammatory cytokines relative to pro-inflammatory cytokines
    S_n = rates[23]                 # strength of PAMPs & DAMPs relative to pro-inflammatory cytokines
    
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
        elif (P + S_n * N - S_a * A) >= P_crit:
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
        elif ((P + S_n * N - S_a * A) > P_crit):
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

        # if count < !NUM_HERE!:
        
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
        # functions should, in theory, never go negative under expected circumstances
    #------------- 3. Output --------------
    
    output = np.array([H_output, N_output, S_output, Q_output, U_output, P_output, A_output, T_output, I_output])
    
    return output


def linear_sim_smooth(parameter_arr, init_value_arr, delta_t, t_final, pathogen_t, pathogen_size):

    '''
    ARGS:
    ---------------
    parameter_arr: Array containing the values for the parameters of the system, array rather than dictionary to maintain compatibility with sensitivity analysis code
    init_value_arr: Array containing the values for the model outputs at time t=0
    delta_t: Timestep size
    t_final: Final timestep to end simulation on
    pathogen_t: Array containing the times at which an external pathogen input will be applied to the system
    pathogen_size: Array containing the size of the external pathogen input to be applied to the system; must match the length of pathogen_t
    
    OUTPUT:
    -----------------
    Array of shape (# of outputs, t_final/delta_t) OR the int '0' if an error occured
    '''

    #--------- 0. Verifying arguments work ---------------
    if len(pathogen_t) != len(pathogen_size):
        return 0
    
    #--------- 1. Initiliazing all values and parameters ----------------

    timesteps = np.arange(0, t_final, delta_t)

    # (X)_t are the variables that will be updated each timestep
    H_t = init_value_arr[0]
    N_t = init_value_arr[1]
    Q_t = init_value_arr[2]
    S_t = init_value_arr[3]
    U_t = init_value_arr[4]
    P_t = init_value_arr[5]
    A_t = init_value_arr[6]

    # total variable parameter count: 25
    g_N = parameter_arr[0]
    k_nq = parameter_arr[1]
    k_ns = parameter_arr[2]
    k_sn = parameter_arr[3]
    S_PQ = parameter_arr[4]
    S_PH = parameter_arr[5]
    S_PS = parameter_arr[6]
    S_AS = parameter_arr[7]
    S_AH = parameter_arr[8]
    S_AU = parameter_arr[9]
    d_s = parameter_arr[10]
    d_p = parameter_arr[11]
    d_a = parameter_arr[12]
    d_q = parameter_arr[13]
    d_u = parameter_arr[14]
    N_half = parameter_arr[15]
    R_MAX = parameter_arr[16]
    R_MIN = parameter_arr[17]
    D_MAX = parameter_arr[18]
    D_MIN = parameter_arr[19]
    R_CRIT = parameter_arr[20]
    D_CRIT = parameter_arr[21]
    S_half = parameter_arr[22]
    U_CRIT = parameter_arr[23]
    Q_CRIT = parameter_arr[24]
    Q_SHIFT = parameter_arr[25]
    U_SHIFT = parameter_arr[26]
    mu_sa_MAX = parameter_arr[27]
    mu_sp_MAX = parameter_arr[28]
    S_N = parameter_arr[29]
    S_A = parameter_arr[30]

    # the following parameters have been fixed for the time being
    N_inf = 2*(10**7)
    
    # the following are the omega functions appearing in the sigmoids, epsilons fixed to 0.01
    omega_R = np.log(((R_MAX-R_MIN)/0.01)-1)/R_CRIT
    omega_D = np.log(((D_MAX-D_MIN)/0.01)-1)/D_CRIT
    omega_U = -1 * np.log(((mu_sa_MAX)/0.01)-1)/(U_CRIT)
    omega_Q = -1 * np.log(((mu_sp_MAX)/0.01)-1)/(Q_CRIT)

    H_output = [H_t]
    N_output = [N_t]
    Q_output = [Q_t]
    S_output = [S_t]
    U_output = [U_t]
    P_output = [P_t]
    A_output = [A_t]
    T_output = [Q_t + S_t + U_t]        # T = combined total of leukocytes (Q + S + U)
    I_output = [P_t + N_t - A_t]        # I = total net inflammation

    count = 1

    #-------------- 2. Run the simulation ---------------------

    H = H_output[0]         # initializing the variables that will be updated each timestep
    N = N_output[0]
    Q = Q_output[0]
    S = S_output[0]
    U = U_output[0]
    P = P_output[0]
    A = A_output[0]
    I = I_output[0]

    for t in timesteps[1:]:
        
        count += 1
        #--------- Mechanism for introducing pathogen input ---------

        for i in range(len(pathogen_t)):
            if count == pathogen_t[i]/delta_t:
                N += pathogen_size[i]

        #------------- 2a. Calculate derivatives -------------------

        R_t = R_MAX - (R_MAX - R_MIN)/(1 + np.exp(omega_R*I))
        D_t = D_MAX - (D_MAX - D_MIN)/(1 + np.exp(omega_D*I))
        mu_sa_t = mu_sa_MAX * 1/(1 + np.exp(omega_U*(I - U_SHIFT)))
        mu_sp_t = mu_sp_MAX * 1/(1 + np.exp(omega_Q*(I - Q_SHIFT)))

        dH = (R_t - D_t)*H
        dN = g_N*N*(1-N/N_inf) - (k_nq*Q + k_ns*S)*(N/(N_half + N))
        dS = D_t*H - (k_sn*N*(S/(S_half+S))) - d_s*S - (mu_sa_t + mu_sp_t)*S
        dQ = mu_sp_t*S - d_q*Q
        dU = mu_sa_t*S - d_u*U
        dP = S_PH*H + S_PS*S + S_PQ*Q - d_p*P
        dA = S_AH*H + S_AS*S + S_AU*U -d_a*A

        # ----------- 2b. Diagnostics ----------------------

        '''if pathogen_t[i]/delta_t - 20 < count < pathogen_t[i]/delta_t + 20:
            print(f'T={count}')
            print(f'mu_sp_t: {mu_sp_t}')
            print(f'mu_sa_t: {mu_sa_t}')'''

        #------------ 2c. Update state variables using linear approximation ------------

        H = H + dH*delta_t
        if H < 0:
            H = 0

        N = N + dN*delta_t
        if N < 0:
            N = 0

        S = S + dS*delta_t
        if S < 0:
            S = 0

        Q = Q + dQ*delta_t
        if Q < 0:
            Q = 0

        U = U + dU*delta_t
        if U < 0:
            U = 0

        P = P + dP*delta_t
        if P < 0:
            P = 0

        A = A + dA*delta_t
        if A < 0:
            A = 0
        
        I = P + S_N*N - S_A*A
        T = S + Q + U

        #-------------- 2d. Append to list -------------------------
        
        H_output.append(H)
        N_output.append(N)
        S_output.append(S)
        Q_output.append(Q)
        U_output.append(U)
        P_output.append(P)
        A_output.append(A)
        T_output.append(T)
        I_output.append(I)

    #----------------3. Return output arrays ---------------------
    # returned output is a matrix of shape (9, t_final/delta_t)
    
    output = [H_output, N_output, S_output, Q_output, U_output, P_output, A_output, T_output, I_output]
    return output


def merge_figures_grid(nRow, nCol, img_width, img_height, exp_num, order, o_names, p_names):
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

    ====================
    Outputs:
    No output; Saves resulting image file, split into an nRow x nCol grid pattern composed of each individual image, into the
    same file location as the source code that called this function
    '''
    path_titles = p_names

    relative_path_arr = np.empty(len(o_names), dtype=f'<U256')    # will contain paths to the images in the same order that the output names appear in o_names arg

    script_dir = os.path.dirname(os.path.abspath(__file__))

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
            print('debug check')
            os.makedirs(filepath)

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

    WARNING: order='second' not fully implemented
    ==========================
    No output, generates time-series figures from given .csv files and stores them in predefined folder locations
    '''
    script_dir = ''
    if filepath=='default':

        script_dir = os.path.dirname(os.path.abspath(__file__))     # file path is relative to source code path

    else:
        script_dir = filepath


    nParam = len(p_names)

    # ========== 1: Load data from csv's into useable format ==========

    master_df = np.zeros((len(o_names), nTimesteps, nParam, nDatapoints))
    # master_df: output -> SIs for all params sorted by time -> SIs sorted by param e.g. to access the sensitivity index of paramater 'z' in the timestep 'y' for output 'x', the index would be df[x, y, z, 1]

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


def calculate_derivatives(values):

    # this function is mostly a copy-paste of the linear_sim function; it is used to return the derivative values for given input values (passed in as the arg 'values')

    g_N = values[0]                  # pathogen growth coefficient per hour
    k_nq = values[1]                 # kill rate of pthogens by active WBCs
    k_ns = values[2]                 # kill rate of pthogens by stable WBCs
    k_tn = values[3]                 # kill rate of stable WBCs by stable leukocytes
    w = values[4]                    # how quickly does E(t) -> 0 (as a function of inflammation I(t)); considering fixing
    P_crit = values[5]               # value of inflammation at which 
    S_a = 1                         # fixed (SA)
    S_n = 1                         # fixed (SA)
    N_inf = 2*10**7                 # fixed (SA)
    S_PQ = values[6]                 # secretion rate of pro by active WBCs
    S_PH = values[7]                 # secretion rate of pro by HSPCs
    S_PS = values[8]                 # secretion rate of pro by stable WBCs
    S_AS = values[9]                 # secretion rate of anti by stable WBCs
    S_AH = values[10]                # secretion rate of anti by HSPCs
    S_AU = values[11]                # secretion rate of anti by immuno-suppressive WBCs
    Immune_start = values[12]        # value of inflammation at which stable->immuno-suppresive WBC conversion begins
    Active_start = values[13]        # value of inflammation at which stable->active WBC conversion begins
    Immune_crit = values[14]         # value of inflammation at which stable->immuno-suppresive WBC reaches its max value
    Active_crit = values[15]         # value of inflammation at which stable->active WBC conversion reaches its max value
    y = values[16]                   # modulates the pathogen killing rate of stable WBCs by the density of stable WBCs
    d_s = values[17]                 # natural decay rate of stable WBCs
    d_p = values[18]                 # natural decay rate of pro-inflammatory cytokines
    d_a = values[19]                 # natural decay rate of anti-inflammatory cytokines
    d_q = values[20]                 # natural decay rate of active WBCs
    d_u = values[21]                 # natural decay rate of immuno_suppressive WBCs
    H = values[22]
    N = values[22]
    P = values[22]
    A = values[22]
    Q = values[22]
    S = values[22]
    U = values[22]
    H_init = values[22]

    "Function for E (exhaustion - dH/dt)"
    
    if (P + S_n * N - S_a * A) < P_crit:
        E = 1
    elif (P + S_n * N - S_a * A) >= P_crit:
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
    elif ((P + S_n * N - S_a * A) > P_crit):
        D = 0.35*H
        
    '''"Determining E_star for H_stable"
    
    if P + S_n*N - S_a*A <= 1.5*P_crit:
        dE_star = 0
        
    else:
        dE_star = -0.0005*(P + S_n*N - S_a*A - 1.5*P_crit)'''
        
    
    
    "BM Niche Renewal contribution (dH/dt)"     #E_star=1
    
    if 1 * H_init > H:
        H_stable = 0.05 * H * (1-(H/(1 * H_init)))
    
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
    

    return [dHdt, dNdt, dPdt, dAdt, dSdt, dQdt, dUdt]
