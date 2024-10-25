import numpy as np


def linear_like(x, C, e):
    '''
    This function acts like f(x)=x for x >> e, approaches C for x << e
    Used in model where ratios are concerned and the denominator can equal 0 which would otherwise create an error

    '''

    return C + (x-C)*(x**2/(x**2 + e))

def beta_model(t, y, parameters, return_terms=False):
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

    # ----------- 3. Calculate functions in derivatives --------------------
    I_t = (P_t*((N_t**k)/(theta_N**k + N_t**k) + 0.25)*((0.5 * (K_t)/(theta_K + K_t) + 1))) / ((A_t*(theta_N**k/(theta_N**k + N_t**k) + 0.25)*(0.75 * (K_t)/(theta_K + K_t) + 1)) + (P_t*((N_t**k)/(theta_N**k + N_t**k) + 0.25)*((0.5 * (K_t)/(theta_K + K_t) + 1))))
    #I_t = P_t / (P_t + A_t)

    #D_I = (HM_t*(1/5)*P_t)/(HM_t + (1/5)*P_t) * I_t              # proportion of proliferating HSPCs differentiating (either symmetric or asymmetric)
    D_I = (HM_t*(1/5)*P_t)/((1/5)*HM_t + P_t) * I_t

    beta = I_t/((I_crit + 10) + I_t)                              # proportion of differentiating proliferating HSPCs asymmetrically differentiating (1 parent HSPC -> 2 daughter WBCs)

    #eta_Q = (1/5) * ((HQ_t * (P_t)) / ((1/5)*HQ_t + P_t)) * (I_t/(I_crit + I_t))
    eta_Q = (1/5) * ((HQ_t * (P_t)) / ((1/5)*HQ_t + P_t)) * I_t

    #eta_M = (1/5) * ((HM_t * A_t) / (1/5*HM_t + A_t)) * ((1/I_t)/(A_crit + (1/I_t)))
    eta_M = (1/5) * ((HM_t * A_t) / (1/5*HM_t + A_t)) * (1 - I_t)

    #print(f'eta_M: {eta_M}')
    # IMPORTANT: These next two functions control how the stable leukocytes compartment (S) upregulate the immuno-suppressive and active compartments (U, Q respectively); relates to our research question

    D_Q = (1/3)*tau_Q * (I_t)/(I_crit + I_t)
    D_U = (1/3)*tau_U * (1-I_t) / (A_crit + (1 - I_t))

    # ----------- debug ------------------



    # ----------- 4. Calculate derivatives --------------------

    dHM_dt = eta_Q - D_I*beta - dH*HM_t - eta_M
    #dHM_dt = eta_Q - D_I*beta - dH*HM_t - eta_M - (0.001*N_t*((1 - (1/5 + dH + 1/5))*HM_t/(0.001*N_t+(1 - (1/5 + dH + 1/5))*HM_t))) left off here; producing "invalid value encountered in scalar divide"

    dHQ_dt = (eta_M - eta_Q) + (2*HQ_t*(1 - (k_H*HQ_t)/linear_like(SCSF_t, 0.1, 0.00001)))
    #dHQ_dt = eta_M - eta_Q + (2*HQ_t*(1 - (k_H*HQ_t)/linear_like(SCSF_t, 0.1, 0.00001))) - (0.000001*N_t*((4/5)*HM_t/(0.000001*N_t+(4/5)*HM_t))) left off here; producing "invalid value encountered in scalar divide"

    #dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t - (k_sn*N_t*(S_t/(k_sn*N_t+S_t)))
    #dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t
    dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t - (k_sn*N_t*((1 - (D_Q + D_U + d_S))*S_t/(k_sn*N_t+(1 - (D_Q + D_U + d_S))*S_t)))   # best so far
    #dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - D_Q*S_t - D_U*S_t - d_S*S_t + (k_sn*N_t*((1 - D_Q + D_U + d_S)*S_t/(k_sn*N_t+(1 - D_Q + D_U + d_S)*S_t)))


    dQ_dt = D_Q*S_t - d_Q*(1 - 0.5*I_t/(2 + I_t))*Q_t

    dU_dt = D_U*S_t - d_U*U_t

    dSCSF_dt = S_SCSF*(K_crit / (K_crit + K_t)) - d_SCSF*SCSF_t

    dP_dt = (S_PH*HM_t + S_PQ*Q_t)*(0.5*I_t/(2 + I_t) + 0.5) + S_PS*S_t - d_P*P_t

    dA_dt = (S_AH*HM_t + S_AU*U_t)*(0.5*I_t/(2 + I_t) + 0.5)  + S_AS*S_t - d_A*A_t

    dK_dt = S_KD*(k_sn*N_t*(S_t/(k_sn*N_t+S_t))) - R_KU*U_t*(K_t/(K_crit + K_t))

    dN_dt = g_N*N_t*(1-(N_t/N_oo)) - (k_nq*Q_t + k_ns*S_t)*(N_t/(N_half + N_t))

    # ---------- 5. Return derivative -------------------

    if return_terms:

        debug_info = np.array([
            eta_Q,
            eta_M, 
            D_I, 
            beta, 
            dH*HM_t,
            D_Q, 
            D_U, 
            (k_sn*N_t*(S_t/(k_sn*N_t+S_t))), 
            g_N*N_t*(1-N_t/N_oo),
            (k_nq*Q_t + k_ns*S_t)*(N_t/(N_half + N_t))
        ])

        return (np.array([dHQ_dt, dHM_dt, dN_dt, dP_dt, dA_dt, dSCSF_dt, dK_dt, dQ_dt, dS_dt, dU_dt]), debug_info)
        
    else:

        return np.array([dHQ_dt, dHM_dt, dN_dt, dP_dt, dA_dt, dSCSF_dt, dK_dt, dQ_dt, dS_dt, dU_dt])
