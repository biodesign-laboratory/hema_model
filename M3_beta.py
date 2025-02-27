import numpy as np


def linear_like(x, C, e):
    '''
    This function acts like f(x)=x for x >> e, approaches C for x << e
    Used in model where ratios are concerned and the denominator can equal 0 which would otherwise create an "invalid value encountered in scalar divide" warning

    '''

    return C + (x-C)*(x**2/(x**2 + e))

def beta_model_3(t, y, parameters):

    '''
    This is model w/ MDSCs
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
    d_MF = parameters['d_MF']              # base decay rate of molecular factorst
    d_S = parameters['d_S']                # base decay rate of S cells
    d_Q = parameters['d_Q']                # base decay rate of Q cells
    d_U = parameters['d_U']                # base decay rate of U cells
    d_P = parameters['d_P']                # base decay rate of pro-inflammatory cytokines
    d_A = parameters['d_A']                # base decay rate of anti-inflammatory cytokines
    d_M = parameters['d_M']                # base decay rate of MDSCs

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
    S_AM = parameters['S_AM']

    S_SCSF = parameters['S_SCSF']              # secretion of stem cell supporting factors (SCSF) by the BM Niche; this is taken to be the maximum rate
    S_MF = parameters['S_MF']                  # upregulation of molecular factors (MF) from surrounding healthy tissue; this is taken to be the maximum rate

    S_KD = parameters['S_KD']                  # Rate of upregulation of tissue damage (in the form of DAMPs) by per cell death
    S_KMD = parameters['S_KMD']                # Rate of upregulation of tissue damaging molecules (e.g. ROS, NO, etc.) by MDSC 
    S_KQ = parameters['S_KQ']                  # Rate of upregulation of tissue damaging molecules (e.g. ROS, NO, etc.) by Q cell 

    # kill rate parameters
    k_sn = parameters['k_sn']                  # N kills S
    k_nq = parameters['k_nq']                  # Q kills N (generally high)
    k_ns = parameters['k_ns']                  # S kills N (generally low)
    k_nm = parameters['k_nm']                  # MDSC kills N (through ROS, generally low)

    # misc parameters
    R_KU = parameters['R_KU']                  # rate at which immuno-suppressive cells heal tissue damage
    I_crit = parameters['I_crit']
    K_crit = parameters['K_crit']              # concentration of DAMPs (K(t)) needed to dampen SCSF production by 0.5x
    k = parameters['k']                        # hill-type coefficient

    psi = parameters['psi']                    # parameter that determines the upper limit of proportion of S cells that can be lost in any one moment
    omega = parameters['omega']                # ratio  (S+P -> MDSC) / (S+P -> Q), see implementation in derivatives below
    alpha = parameters['alpha']                # Part of determining maximum velocity of HQ + P -> HM and HM + P -> S, see implementation in eta_Q

    beta_N = parameters['beta_N']              # how quickly N upregulates K

    C_QM = parameters['C_QM']                  # consumption rate of molecular factors by pro-inflammatory WBCs 
    C_MDM = parameters['C_MDM']                # consumption rate of molecular factors by MDSCs
    C_UM = parameters['C_UM']                  # consumption rate of molecular factors by anti-inflammatory WBCs
    C_UP = parameters['C_UP']                  # consumption rate of pro-inflammatory cytokines by anti-inflammatory WBCs



    # ----------- 2. Load variable values --------------------

    HQ_t = y[0]             # Quiescent HSPCs
    HM_t = y[1]             # Mobilized HSPCs
    N_t = y[2]              # Pathogen (PAMPs)
    P_t = y[3]              # Pro-inflammatory cytokines
    A_t = y[4]              # Anti-inflammatory cytokines
    SCSF_t = y[5]           # Stem Cell Supporting Factors (secreted by BM Niche tissue)
    K_t = y[6]              # Tissue Damage (DAMPs)
    Q_t = y[7]              # Active Leukocytes (e.g. M1 Macrophages, NK Cells, Neutrophils, etc.)
    S_t = y[8]              # Pluri-potent, immature WBCs ready for differentiation
    U_t = y[9]              # Immuno-suppressive Leukocytes (e.g. M2 Macrophages, T-cells, etc.)
    MDSC_t = y[10]          # MDSCs
    MF_t = y[11]            # Molecular Factors (e.g. L-arginine, iron, etc.)

    # ----------- 3a. Calculate functions in derivatives --------------------
    amp_P_by_N = (N_t**k)/(theta_N**k + N_t**k) + 0.25      # amplifies pro-inflammatory signals by pathogen concentration
    amp_P_by_K = 0.5 * (K_t)/(theta_K + K_t) + 1            # amplifies pro-inflammatory signals by tissue damage
    # damp_A_by_N = theta_N**k/(theta_N**k + N_t**k) + 0.25   # dampens the anti-inflammatory signals by pathogen concentration
    amp_A_by_K = 0.75 * (K_t)/(theta_K + K_t) + 1           # amplifies the anti-inflammatory signals by tissue damage

    amp_A_by_SCSF = 0.75 * (SCSF_t)/(theta_K + SCSF_t) + 1

    I_t = (P_t * amp_P_by_N * amp_P_by_K) / ((A_t * amp_A_by_K) + (P_t * amp_P_by_N * amp_P_by_K))

    #I_H = (P_t * amp_P_by_N * amp_P_by_K) / ((A_t * amp_A_by_K) + (P_t * amp_P_by_N * amp_P_by_K) + SCSF_t)
    I_H = (P_t * amp_P_by_N * amp_P_by_K) / ((A_t * amp_A_by_K + amp_A_by_SCSF) + (P_t * amp_P_by_N * amp_P_by_K))

    D_I = (HM_t*(alpha*I_t)*P_t)/((alpha*I_t)*HM_t + P_t)      # proportion of proliferating HSPCs differentiating (either symmetric or asymmetric)

    beta = I_t/(I_crit + I_t)                              # proportion of differentiating proliferating HSPCs asymmetrically differentiating (1 parent HSPC -> 2 daughter WBCs)

    eta_Q = (alpha*I_H) * ((P_t * HQ_t) / ((alpha*I_H)*HQ_t + P_t))    # currently used term

    downregulate_S = (tau_Q*P_t + tau_U*A_t + k_sn*N_t)*psi*S_t / (tau_Q*P_t + tau_U*A_t + k_sn*N_t + psi*S_t)      # controls number of S cells lost due to interaction with P, A, or N
    I_S = tau_Q*P_t*amp_P_by_N*amp_P_by_K + tau_U*A_t*amp_A_by_K + k_sn*N_t
    # I_S = tau_Q*P_t + tau_U*A_t + k_sn*N_t          # w/out amp functions      


    # ----------- 3b. debug ------------------  


    # ----------- 4. Calculate derivatives --------------------

    dHM_dt = eta_Q - D_I*beta - dH*HM_t

    dHQ_dt = (2*HQ_t*(1 - (k_H*HQ_t)/linear_like(SCSF_t, 0.1, 0.00001))) - eta_Q

    dS_dt = (D_I*(1 - beta) + 2*D_I*beta) - downregulate_S

    dQ_dt = (omega) * downregulate_S * (tau_Q*P_t*amp_P_by_N*amp_P_by_K) / (I_S) * ((1/2*MF_t) / (1/2*MF_t + C_QM*Q_t + C_MDM*MDSC_t + C_UM*U_t)) - d_Q*(1 - 0.5*I_t/(0.7 + I_t))*Q_t      # Why multiply by 1/2? Because S + P -> Q AND MDSC, so I simplify here and assume half become Q, half become MDSCs
    # dQ_dt = (3/4) * downregulate_S * (tau_Q*P_t) / (I_S) * ((1/2*MF_t) / (1/2*MF_t + C_QM*Q_t + C_MDM*MDSC_t + C_UM*U_t)) - d_Q*(1 - 0.5*I_t/(0.7 + I_t))*Q_t           # w/out amp functions

    dU_dt = downregulate_S * (tau_U*A_t*amp_A_by_K) / (I_S) - d_U*U_t
    # dU_dt = downregulate_S * (tau_U*A_t) / (I_S) - d_U*U_t           # w/out amp functions

    dSCSF_dt = S_SCSF*(K_crit / (K_crit + K_t)) - d_SCSF*SCSF_t

    dP_dt = (S_PH*HM_t + S_PQ*Q_t)*(0.8*I_t + 0.2) + S_PS*S_t - d_P*P_t - (C_UP*U_t * (1-d_P)*P_t)/(C_UP*U_t + (1-d_P)*P_t)

    dA_dt = S_AM*MDSC_t + S_AH*HM_t + S_AU*U_t  + S_AS*S_t - d_A*A_t
    
    # dK_dt = S_KD*((k_sn*N_t*(S_t/(k_sn*N_t+S_t)))*(k_sn*N_t / I_S)) + S_KMD*MDSC_t*((P_t + A_t + N_t)/(1/2*MDSC_t + P_t + A_t + N_t)) + S_KQ*Q_t - (R_KU*U_t*K_t / (R_KU*U_t + K_t))
    dK_dt = beta_N*N_t + S_KD*((k_sn*N_t*(S_t/(k_sn*N_t+S_t)))*(k_sn*N_t / I_S)) + S_KMD*MDSC_t*((P_t + A_t + N_t)/(1/2*MDSC_t + P_t + A_t + N_t)) + S_KQ*Q_t - (R_KU*U_t*K_t / (R_KU*U_t + K_t))
    # dK_dt = beta_N*N_t + S_KD*((k_sn*N_t*(S_t/(k_sn*N_t+S_t)))*(k_sn*N_t / I_S) + eta_Q*(beta_N*N_t / (beta_N*N_t + P_t))) + S_KMD*MDSC_t*((P_t + A_t + N_t)/(1/2*MDSC_t + P_t + A_t + N_t)) + S_KQ*Q_t - (R_KU*U_t*K_t / (R_KU*U_t + K_t))

    dN_dt = g_N*N_t*(1-(N_t/N_oo)) - (k_nq*Q_t + k_ns*S_t + k_nm*MDSC_t )*(N_t/(N_half + N_t))

    dMDSC_dt = (1-omega) * downregulate_S * (tau_Q*P_t*amp_P_by_N*amp_P_by_K) / linear_like(I_S, 0.1, 0.00001) - d_M*( 0.5*Q_t**k / linear_like(1/2*MDSC_t**k + Q_t**k, 0.1, 0.00001) + 0.5)*MDSC_t
    # dMDSC_dt = (1/4) * downregulate_S * (tau_Q*P_t) / (I_S) - d_M*( 0.5*Q_t**k / (1/2*MDSC_t**k + Q_t**k) + 0.5)*MDSC_t

    # dTDM_dt = S_TM/3 * ( ((MDSC_t * eps_3 * P_t) / (MDSC_t + eps_3*P_t)) + ((MDSC_t * eps_4 * A_t) / (MDSC_t + eps_4*A_t)) + ((MDSC_t * N_t) / (MDSC_t + N_t))) - d_TDM*TDM_t     Unnecessary ?

    dMF_t = S_MF*(K_crit/(K_crit + K_t)) - d_MF*MF_t

    # ---------- 5. Return derivative -------------------

    return np.array([dHQ_dt, dHM_dt, dN_dt, dP_dt, dA_dt, dSCSF_dt, dK_dt, dQ_dt, dS_dt, dU_dt, dMDSC_dt, dMF_t])

