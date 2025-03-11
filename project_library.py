from scipy.integrate import solve_ivp
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
    R_SHIFT = parameter_arr[31]
    D_SHIFT = parameter_arr[32]

    # the following parameters have been fixed for the time being
    N_inf = 2*(10**7)
    
    # the following are the omega functions appearing in the sigmoids, epsilons fixed to 0.01
    omega_R = np.log(((R_MAX-R_MIN)/0.01)-1)/R_CRIT
    omega_D = np.log(((D_MAX-D_MIN)/0.01)-1)/D_CRIT
    omega_U = np.log(((mu_sa_MAX)/0.01)-1)/(U_CRIT)
    omega_Q = np.log(((mu_sp_MAX)/0.01)-1)/(Q_CRIT)

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

        R_t = R_MAX - (R_MAX - R_MIN)/(1 + np.exp(omega_R*(I - R_SHIFT)))
        D_t = D_MAX - (D_MAX - D_MIN)/(1 + np.exp(omega_D*(I - D_SHIFT)))
        mu_sa_t = mu_sa_MAX * 1/(1 + np.exp(omega_U*(I - U_SHIFT)))
        mu_sp_t = mu_sp_MAX * 1/(1 + np.exp(omega_Q*(I - Q_SHIFT)))

        #print(D_t*H)

        dH = (R_t - D_t)*H
        dN = g_N*N*(1-N/N_inf) - (k_nq*Q + k_ns*S)*(N/(N_half + N))
        dS = D_t*H - (k_sn*(10**4 * N / (N_half + N))*(S/(S_half+S))) - d_s*S - (mu_sa_t + mu_sp_t)*S
        dQ = mu_sp_t*S - d_q*Q
        dU = mu_sa_t*S - d_u*U
        dP = S_PH*H + S_PS*S + S_PQ*Q - d_p*P
        dA = S_AH*H + S_AS*S + S_AU*U -d_a*A

        # ----------- 2b. Diagnostics ----------------------

        '''if pathogen_t[i]/delta_t - 20 < count < pathogen_t[i]/delta_t + 20:
            print(f'T={count}')
            print(f'mu_sp_t: {mu_sp_t}')
            print(f'mu_sa_t: {mu_sa_t}')'''
        '''if count < 12000 and count > 9500 and count % 20 == 0:
            print('===================================')
            print(f"Timestep {t} : R(t)={R_t}")
            print(f"Timestep {t} : D(t)={D_t}")'''

        #------------ 2c. Update state variables using linear approximation ------------

        if H <= 100:
            H = 0
            #print(f'{t}:H <= 100')
        else:
            H = H + dH*delta_t

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
        #print(f"I(t)={I}")
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


def calculate_derivatives_hybrid(values):

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


def calculate_derivatives_smooth(t, y, parameters):

    '''
    This function is mostly a copy-paste of the derivate section of linear_sim_smooth made compatible with solve_ivp
    ARGS:
    t: scalar (needed by solve_ivp)
    y: ndarray with the form (H(t), N(t), P(t), A(t), Q(t), S(t), U(t))
    parameters: dictionary of parameter values
    ==============================
    Outputs derivative values given inputs
    '''
    # -------- 1. Unpack dictionary, define relevant functions ----------------------------
    g_n = parameters['g_n']
    k_nq = parameters['k_nq']
    k_ns = parameters['k_ns']
    k_sn = parameters['k_sn']
    S_PQ = parameters['S_PQ']
    S_PH = parameters['S_PH']
    S_PS = parameters['S_PS']
    S_AS = parameters['S_AS']
    S_AH = parameters['S_AH']
    S_AU = parameters['S_AU']
    d_s = parameters['d_s']
    d_p = parameters['d_p']
    d_a = parameters['d_a']
    d_q = parameters['d_q']
    d_u = parameters['d_u']
    N_half = parameters['N_half']
    R_MAX = parameters['R_MAX']
    R_MIN = parameters['R_MIN']
    D_MAX = parameters['D_MAX']
    D_MIN = parameters['D_MIN']
    R_CRIT = parameters['R_CRIT']
    D_CRIT = parameters['D_CRIT']
    S_half = parameters['S_half']
    U_CRIT = parameters['U_CRIT']
    Q_CRIT = parameters['Q_CRIT']
    Q_SHIFT = parameters['Q_SHIFT']
    U_SHIFT = parameters['U_SHIFT']
    mu_SA_MAX = parameters['mu_SA_MAX']
    mu_SP_MAX = parameters['mu_SP_MAX']
    S_N = parameters['S_N']
    S_A = parameters['S_A']
    R_SHIFT = parameters['R_SHIFT']
    D_SHIFT = parameters['D_SHIFT']

    N_inf = 2 * (10**7)     # 20 million, fixed

    omega_R = np.log(((R_MAX-R_MIN)/0.01)-1)/R_CRIT
    omega_D = np.log(((D_MAX-D_MIN)/0.01)-1)/D_CRIT
    omega_U = np.log(((mu_SA_MAX)/0.01)-1)/(U_CRIT)
    omega_Q = np.log(((mu_SP_MAX)/0.01)-1)/(Q_CRIT)

    # ---------- 2. Calculate derivatives -------------------
    H_t = y[0]
    N_t = y[1]
    P_t = y[2]
    A_t = y[3]
    Q_t = y[4]
    S_t = y[5]
    U_t = y[6]

    I_t = P_t + S_N*N_t - S_A*A_t
    
    R_t = R_MAX - (R_MAX - R_MIN)/(1 + np.exp(omega_R*(I_t - R_SHIFT)))
    D_t = D_MAX - (D_MAX - D_MIN)/(1 + np.exp(omega_D*(I_t - D_SHIFT)))
    mu_sa_t = mu_SA_MAX * 1/(1 + np.exp(omega_U*(I_t - U_SHIFT)))
    mu_sp_t = mu_SP_MAX * 1/(1 + np.exp(omega_Q*(I_t - Q_SHIFT)))

    dH = (R_t - D_t)*H_t
    dN = g_n*N_t*(1-N_t/N_inf) - (k_nq*Q_t + k_ns*S_t)*(N_t/(N_half + N_t))
    dS = D_t*H_t - (k_sn*N_t*(S_t/(N_t+S_t))) - d_s*S_t - (mu_sa_t + mu_sp_t)*S_t
    dQ = mu_sp_t*S_t - d_q*Q_t
    dU = mu_sa_t*S_t - d_u*U_t
    dP = S_PH*H_t + S_PS*S_t + S_PQ*Q_t - d_p*P_t
    dA = S_AH*H_t + S_AS*S_t + S_AU*U_t -d_a*A_t

    return [dH, dN, dP, dA, dQ, dS, dU]


def calculate_derivatives_2(t, y, parameters):
    '''
    t: scalar (needed by solve_ivp)
    y: ndarray with the form (H(t), N(t), S(t), Q(t), U(t), P(t), A(t), K(t))
    parameters: dictionary of parameter values
    ==============================
    Outputs derivative values given inputs
    This is the calculate_derivatives function for model version #2
    '''

    # -------- 1. Unpack dictionary, define relevant functions ----------------------------
    g_n = parameters['g_n']
    k_nq = parameters['k_nq']
    k_ns = parameters['k_ns']
    k_sn = parameters['k_sn']
    S_PQ = parameters['S_PQ']
    S_PH = parameters['S_PH']
    S_PS = parameters['S_PS']
    S_AS = parameters['S_AS']
    S_AH = parameters['S_AH']
    S_AU = parameters['S_AU']
    d_s = parameters['d_s']
    d_p = parameters['d_p']
    d_a = parameters['d_a']
    d_q = parameters['d_q']
    d_u = parameters['d_u']
    R_MAX = parameters['R_MAX']
    D_MAX = parameters['D_MAX']
    S_N = parameters['S_N']
    R_half = parameters['R_half']
    D_half = parameters['D_half']
    S_K = parameters['S_K']
    I_half = parameters['I_half']
    U_half = parameters['U_half']
    Q_half = parameters['Q_half']
    U_MAX = parameters['U_MAX']
    Q_MAX = parameters['Q_MAX']
    d_k = parameters['d_k']
    K_half = parameters['K_half']
    S_KS = parameters['S_KS']
    tau = parameters['tau']
    tau_half = parameters['tau_half']
    H_oo = parameters['H_oo']

    N_inf = 2 * (10**7)     # 20 million, fixed
    N_half = 1000           # fixed
    S_half = 1000           # fixed

    # ------------- 2. Set variable values -----------------
    H_t = y[0]
    N_t = y[1]
    S_t = y[2]
    Q_t = y[3]
    U_t = y[4]
    P_t = y[5]
    A_t = y[6]
    K_t = y[7]

    f_A = I_half / (I_half + A_t)
    I_t = f_A * (P_t + S_N*N_t + S_K*K_t)

    # -------------- 3. Calculate derivatives ---------------
    R_i = R_MAX * N_t / (R_half + N_t) + 0.1                  # self-renewal (dh/dt)
    D_i = D_MAX * N_t / (D_half + N_t) + 0.1                  # differentiation (dh/dt)

    S_PH_i = tau * (I_t**4 / (tau_half**4 + I_t**4))      # HSPC pro-inflammatory cytokine secretion modulated by inflammation, fixed hill coefficient to k=4

    dH_dt = (R_i - D_i)*H_t + ((1-(H_t / H_oo))*H_t)*(I_t/(I_t + 100)) - 0.05*H_t
    dN_dt = g_n*N_t*(1-(N_t/N_inf)) - (k_nq*Q_t + k_ns*S_t)*(N_t/(N_half + N_t))
    dS_dt = D_i*H_t - k_sn*((10**4 * N_t)/(N_half + N_t))*(S_t/(S_half + S_t)) - U_MAX*(I_t/(U_half + I_t))*S_t - Q_MAX*(I_t/(Q_half + I_t))*S_t - d_s*S_t
    dQ_dt = Q_MAX*(I_t/(Q_half + I_t))*S_t - d_q*Q_t
    dU_dt = U_MAX*(I_t/(U_half + I_t))*S_t - d_u*U_t
    dP_dt = S_PH_i*H_t + S_PS*S_t + S_PQ*Q_t - d_p*P_t
    dA_dt = S_AH*H_t + S_AS*S_t + S_AU*U_t - d_a*A_t
    dK_dt = S_KS*(k_sn*((10**4 * N_t)/(N_half + N_t))*(S_t/(S_half + S_t)) + ((10**4 * S_t)/(S_half + S_t))*(I_t/(K_half + I_t))) - d_k*K_t

    # ------------- 4. Diagnostics -----------------------------
    print(f'Timestep {t}: {k_sn*((10**4 * N_t)/(N_half + N_t))*(S_t/(S_half + S_t))}')

    # ------------- 5. Return derivatives ----------------------

    return [dH_dt, dN_dt, dS_dt, dQ_dt, dU_dt, dP_dt, dA_dt, dK_dt]


def linear_like(x, C, e):
    '''
    This function acts like f(x)=x for x >> e, approaches C for x << e
    Used in model where ratios are concerned and the denominator can equal 0 which would otherwise create an error

    '''

    return C + (x-C)*(x**2/(x**2 + e))


def model_2_derivatives(t, y, parameters):
    '''
    ARGS:
    - parameters : dictionary of relevant parameter values
    - y : model variables
    - t : to make compatible with solve_ivp
    OUTPUTS:
    - derivatives : array of derivative values
    '''
    # ----------- 1. Load parameter values --------------------

    # HSPC parameters
    k_H = parameters['k_H']                 # number of SCSFs consumed per fully renewing proliferating HSPC
    dH = parameters['dH']                   # decay of proliferating HSPCs owing to their sped up cycling

    # sensitivity parameters

    theta_N = parameters['theta_N']
    theta_K = parameters['theta_K']
    tau_Q = parameters['tau_Q']
    tau_U = parameters['tau_U']


    # decay rate parameters

    d_SCSF = parameters['d_SCSF']          # base decay rate of SCSFs
    d_S = parameters['d_S']                # base decay rate of S cells
    d_Q = parameters['d_Q']                # base decay rate of Q cells
    d_U = parameters['d_U']                # base decay rate of U cells
    d_P = parameters['d_P']                # base decay rate of pro-inflammatory cytokines
    d_A = parameters['d_A']                # base decay rate of anti-inflammatory cytokines

    # S_HP = parameters['S_HP']

    # pathogen parameters
    g_N = parameters['g_N']
    N_oo = parameters['N_oo']
    N_half = parameters['N_half']

    # Secretion rate parameters

    S_PH = parameters['S_PH']                  # base rate of P by HM
    S_PS = parameters['S_PS']                  # base rate of P by S
    S_PQ = parameters['S_PQ']                  # base rate of P by Q

    S_AU = parameters['S_AU']
    S_AH = parameters['S_AH']
    S_AS = parameters['S_AS']

    S_SCSF = parameters['S_SCSF']              # secretion of stem cell supporting factors (SCSF) by the BM Niche

    S_KD = parameters['S_KD']                  # amount of DAMPs secreted per cell death

    # kill rate parameters
    k_sn = parameters['k_sn']                  # N kills S
    k_nq = parameters['k_nq']                  # Q kills N (generally high)
    k_ns = parameters['k_ns']                  # S kills N (generally low)

    # misc parameters
    R_KU = parameters['R_KU']                  # rate at which immuno-suppressive cells heal tissue damage
    I_crit = parameters['I_crit']
    A_crit = parameters['A_crit']
    K_crit = parameters['K_crit']              # concentration of DAMPs (K(t)) needed to dampen SCSF production by 0.5x
    k = parameters['k']                        # hill-type coefficient

    psi = parameters['psi']                    # new experimental term


    # ----------- 2. Load variable values --------------------

    HQ_t = y[0]             # Quiescent HSPCs
    HM_t = y[1]             # Mobilized HSPCs
    N_t = y[2]              # Pathogen (PAMPs)
    P_t = y[3]              # Pro-inflammatory cytokines
    A_t = y[4]              # Anti-inflammatory cytokines
    SCSF_t = y[5]           # Stem Cell Supporting Factors (secreted by BM Niche tissue)
    K_t = y[6]              # Tissue Damage (DAMPs)
    Q_t = y[7]              # Active Leukocytes (e.g. M1 Macrophages, NK Cells, Neutrophils, etc.)
    S_t = y[8]              # Stable Leukocytes (i.e. steady-state leukocytes)
    U_t = y[9]              # Immuno-suppressive Leukocytes (e.g. M2 Macrophages, T-cells, etc.)

    # ----------- 3a. Calculate functions in derivatives --------------------
    amp_P_by_N = (N_t**k)/(theta_N**k + N_t**k) + 0.25      # amplifies pro-inflammatory signals by pathogen concentration
    amp_P_by_K = 0.5 * (K_t)/(theta_K + K_t) + 1            # amplifies pro-inflammatory signals by tissue damage
    damp_A_by_N = theta_N**k/(theta_N**k + N_t**k) + 0.25   # dampens the anti-inflammatory signals by pathogen concentration
    amp_A_by_K = 0.75 * (K_t)/(theta_K + K_t) + 1           # amplifies the anti-inflammatory signals by tissue damage

    I_t = (P_t * amp_P_by_N * amp_P_by_K) / ((A_t * damp_A_by_N * amp_A_by_K) + (P_t * amp_P_by_N * amp_P_by_K))

    #D_I = (HM_t*(1/5)*P_t)/(HM_t + (1/5)*P_t) * I_t              # proportion of proliferating HSPCs differentiating (either symmetric or asymmetric)
    #D_I = (HM_t*(1/5)*P_t)/((1/5)*HM_t + P_t) * I_t
    D_I = (HM_t*(1/5)*P_t)/((1/5)*HM_t + P_t)

    beta = I_t/(I_crit + I_t)                              # proportion of differentiating proliferating HSPCs asymmetrically differentiating (1 parent HSPC -> 2 daughter WBCs)

    #eta_Q = (1/5) * ((HQ_t * (P_t)) / ((1/5)*HQ_t + P_t)) * (I_t/(I_crit + I_t))
    eta_Q = (I_t) * ((HQ_t * (P_t)) / ((I_t)*HQ_t + P_t))

    #eta_M = (1/5) * ((HM_t * A_t) / (1/5*HM_t + A_t)) * ((1/I_t)/(A_crit + (1/I_t)))
    eta_M = (1-I_t) * ((HM_t * A_t) / ((1-I_t)*HM_t + A_t))

    #print(f'eta_M: {eta_M}')
    # IMPORTANT: These next two functions control how the stable leukocytes compartment (S) upregulate the immuno-suppressive and active compartments (U, Q respectively); relates to our research question

    '''D_Q = (1/3)*tau_Q * (I_t)/(I_crit + I_t)
    D_U = (1/3)*tau_U * (1-I_t) / (A_crit + (1 - I_t))'''

    # ----------- 3b. debug ------------------
    # new terms to replace D_Q, D_U, and N + S -> __ in dS/dt
    downregulate_S = (tau_Q*P_t + tau_U*A_t + k_sn*N_t)*psi*S_t / (tau_Q*P_t + tau_U*A_t + k_sn*N_t + psi*S_t)
    I_S = tau_Q*P_t*amp_P_by_N*amp_P_by_K + tau_U*A_t*amp_A_by_K*damp_A_by_N + k_sn*N_t


    # ----------- 4. Calculate derivatives --------------------
    # do not delete commented out equations, this is useful for keeping a record of changes in case something breaks

    dHM_dt = eta_Q - D_I*beta - dH*HM_t - eta_M
    #dHM_dt = eta_Q - D_I*beta - dH*HM_t - eta_M - (0.001*N_t*((1 - (1/5 + dH + 1/5))*HM_t/(0.001*N_t+(1 - (1/5 + dH + 1/5))*HM_t))) left off here; producing "invalid value encountered in scalar divide"

    dHQ_dt = (eta_M - eta_Q) + (2*HQ_t*(1 - (k_H*HQ_t)/linear_like(SCSF_t, 0.1, 0.00001)))
    #dHQ_dt = eta_M - eta_Q + (2*HQ_t*(1 - (k_H*HQ_t)/linear_like(SCSF_t, 0.1, 0.00001))) - (0.000001*N_t*((4/5)*HM_t/(0.000001*N_t+(4/5)*HM_t))) left off here; producing "invalid value encountered in scalar divide"

    #dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t - (k_sn*N_t*(S_t/(k_sn*N_t+S_t)))
    #dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t
    #dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t - (k_sn*N_t*((1 - (D_Q + D_U + d_S))*S_t/(k_sn*N_t+(1 - (D_Q + D_U + d_S))*S_t)))   # best so far
    #dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t + (k_sn*N_t*((1 - D_Q + D_U + d_S)*S_t/(k_sn*N_t+(1 - D_Q + D_U + d_S)*S_t)))

    dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - downregulate_S      # new term

    # dQ_dt = D_Q*S_t - d_Q*(1 - 0.5*I_t/(2 + I_t))*Q_t
    dQ_dt = downregulate_S * (tau_Q*P_t*amp_P_by_N*amp_P_by_K) / (I_S) - d_Q*(1 - 0.5*I_t/(2 + I_t))*Q_t      # new term

    # dU_dt = D_U*S_t - d_U*U_t
    dU_dt = downregulate_S * (tau_U*A_t*amp_A_by_K*damp_A_by_N) / (I_S) - d_U*U_t         # new term

    dSCSF_dt = S_SCSF*(K_crit / (K_crit + K_t)) - d_SCSF*SCSF_t

    dP_dt = (S_PH*HM_t + S_PQ*Q_t)*(0.5*I_t/(2 + I_t) + 0.5) + S_PS*S_t - d_P*P_t

    dA_dt = (S_AH*HM_t + S_AU*U_t)*(0.5*I_t/(2 + I_t) + 0.5)  + S_AS*S_t - d_A*A_t

    dK_dt = S_KD*(k_sn*N_t*(S_t/(k_sn*N_t+S_t))) - R_KU*U_t*(K_t/(K_crit + K_t))

    dN_dt = g_N*N_t*(1-(N_t/N_oo)) - (k_nq*Q_t + k_ns*S_t)*(N_t/(N_half + N_t))

    # ---------- 5. Return derivative -------------------

    return np.array([dHQ_dt, dHM_dt, dN_dt, dP_dt, dA_dt, dSCSF_dt, dK_dt, dQ_dt, dS_dt, dU_dt])


def lin_sim(ODE_eq, parameters, init_y, t_final, delta_t, ext_stimuli, ext_stim_method, return_Inflammation = True, return_derivatives=False):
    '''
    ARGS:
    ODE_eq : Function argument, first-order derivatives of system output

    parameters : dictionary of relevant parameters, must be compatible with ODE_eq otherwise an error will most likely occur

    init_y : Initial state of the system, init_y=(y_1(0), y_2(0), ..., y_n(0))

    t_final : Final timestep to reach

    delta_t : Timestep size

    ext_stimuli : np.matrix() -> rows = model output type, columns = specific timestep, array of external inputs to the system; MUST match output matrix in shape

    ext_stim_method: vector of strings, 'ADD' = add external stimuli matrix element wise to model output, 'OVERWRITE' = overwrite calculated model output with ext_stimuli input, this argument
                     must be a vector of strings with the same length as init_y (the number of outputs), each element will denote how to handle the external stimuli input for each matching
                     specific model output (e.g. ext_stim_method[i] == 'ADD' means that for each timestep, the loop will add ext_stimuli[i, timestep] to model_output[i, timestep])

    return_Inflammation : Whether model_output should include the state variable I(t)

    return_derivatives : Whether tuple returned by function should include derivative information, will appear in data[1]

    OUTPUT:
    data : 3 tuple of form (model_output, derivatives_output, debug_output); if return_derivatives = False, derivative_outputs = 0 (same applies to debug_output)
    ''' 
    def calculate_I(P_t, A_t, K_t, N_t, theta_N, theta_K, k):

        return (P_t*(N_t**k/(theta_N**k + N_t**k) + 0.25)*((0.5 * (K_t)/(theta_K + K_t) + 1))) / (A_t*(theta_N**k/(theta_N**k + N_t**k) + 0.25)*(0.75 * (K_t)/(theta_K + K_t) + 1) + P_t*(N_t**k/(theta_N**k + N_t**k) + 0.25)*((0.5 * (K_t)/(theta_K + K_t) + 1)))

    timesteps = np.arange(stop=t_final, step=delta_t)
    model_output = np.zeros((len(init_y), len(timesteps)))
    model_output[:, 0] = init_y
    derivative_output = np.zeros((len(init_y), len(timesteps)))

    if return_Inflammation:
        inflammation_arr = []

    '''if debug_mode:
        debug_output = np.zeros((10, len(timesteps)))
    else:
        debug_output = 0'''

    for i, t in enumerate(timesteps[1:]):
        
        ys = model_output[:, i]

        derivatives = ODE_eq(t, ys, parameters)
        '''if debug_mode:
            derivative_output [:, i] = derivatives[0]
            debug_output = derivatives[1]
            model_output[:, i+1] = ys + delta_t*derivatives[0]        # Euler approximation step here
        else:'''
        derivative_output [:, i] = derivatives
        model_output[:, i+1] = ys + delta_t*derivatives         # Euler approximation step here

        #model_output[:, i+1] = ys + delta_t*derivative_output         # Euler approximation step here

        if return_Inflammation:
            
            if i == 0:
                
                inflammation_arr.append(calculate_I(model_output[3, i], model_output[4, i], model_output[6, i], model_output[2, i], parameters['theta_N'], parameters['theta_K'], parameters['k']))
            
            inflammation_arr.append(calculate_I(model_output[3, i+1], model_output[4, i+1], model_output[6, i+1], model_output[2, i+1], parameters['theta_N'], parameters['theta_K'], parameters['k']))

        for k, method in enumerate(ext_stim_method):
            
            if method == 'ADD':
                model_output[k, i+1] = model_output[k, i+1] + ext_stimuli[k, i+1]
            elif method == 'OVERWRITE':
                model_output[k, i+1] = ext_stimuli[k, i+1]

    if return_Inflammation:
        model_output = np.vstack((model_output, inflammation_arr))

    data = (model_output, derivative_output)

    return data
            


def calculate_I(P_t, A_t, K_t, N_t, theta_N, theta_K, k):
    
    return (P_t*(N_t**k/(theta_N**k + N_t**k) + 0.25)*((0.5 * (K_t)/(theta_K + K_t) + 1))) / (A_t*(theta_N**k/(theta_N**k + N_t**k) + 0.25)*(0.75 * (K_t)/(theta_K + K_t) + 1) + P_t*(N_t**k/(theta_N**k + N_t**k) + 0.25)*((0.5 * (K_t)/(theta_K + K_t) + 1)))


def event_function(t, y, space=None, delta_time=None):
    """
    trigger event when t-delta_time = 0
    space is just a dummy variable. the number of parameters
    after t,y need to match the number of terms in the "args" keyword arg
    in the solve_ivp call.
    """
    
    return t - delta_time

    
def lin_sim_scipy(ODE_eq, parameters, y0, tf, dt, stim_time, stim_size):
    
    '''
    lin_sim for scipy
    
    ARGS:
    ODE_eq : Function argument, first-order derivatives of system output
    parameters : dict, must be compatible with ODE_eq
    init_y : Initial state of the system, init_y=(y_1(0), y_2(0), ..., y_n(0))
    t_final : Final timestep to reach
    delta_t : Timestep size

    stim_time : time(s) of delta-function
    stim_size : stimulus size(s)
    Start with 1
    I will add heaviside inputs next.

    OUTPUT:
    solution
    '''

    event_function.terminal = True

    t = np.arange(0,tf,dt)
    
    out = solve_ivp(ODE_eq,(t[0],t[-1]),y0,
                    args=(parameters,stim_time),
                    events=event_function)

    y2 = out.y[:,-1]
    y2[2] += stim_size
    t2 = np.arange(out.t[-1]+dt,tf,dt)
    print(out.t[-1],tf)

    out2 = solve_ivp(ODE_eq,(t2[0],t2[-1]),y2,
                     args=(parameters,stim_time),
                     method='LSODA')

    data = np.concatenate([out.y.T,out2.y.T])
    t = np.concatenate([out.t,out2.t])
    print(out.t[-1],out2.t[0],out2.t[-1])

    return t,data.T
            
