#Smooth Model Analysis
using ReachabilityAnalysis, ModelingToolkit, Plots, OrdinaryDiffEq, DifferentialEquations

# Extract rates into individual variables
    w = 0.0005
    P_crit = 7000
    A_crit = 1000
    t_half_leukocytes = 7
    t_half_P = 4.1375
    t_half_A = 7
    t_double_P = 15
    gamma = 0.0000001
    S_a = 1
    S_n = 4
    N_inf = 20000000
    S_PQ = 5
    S_PH = 1
    S_PS = 2
    S_AS = 1
    S_AH = 0.5
    S_AU = 7
    theta_ps = 10000000
    theta_ar = 10000000
    theta_AS = 10000000
    theta_UP = 10000000
    Immune_start = 500
    Active_start = 500
    Immune_crit = 2500
    Active_crit = 2500
    y = 2500
    d_s = 0.25
    d_p = 1/4
    d_a = 0.8
    d_q = 0.3
    d_u = 0.6
    g_n = 1.2
    K_PS = 0.7
    K_AQ = 0.7
    K_PU = 0.7
    K_AS = 0.7
    k_nq = 5
    k_ns = 2
    k_sn = 3
    k_tn = 0.33
    H = 1000
    N_func = 0.0
    P = 600
    A = 1160
    T = 7000
    a = 0.00
    b = 1
    e = 0.00
    E_star = 1
    H_0 = 1000
    N_0 = 0.0
    P_0 = 600
    A_0 = 1160
    T_0 = 7000
    a_0 = 0.00
    b_0 = 1
    e_0 = 0.00
    E_star_0 = 1
    Q = a * T
    S = b * T
    U = e * T

#Variables
epsilon_R = 0.001
epsilon_D = 0.001
R_max = 0.35
R_min = 0.1
D_max = 0.45
D_min = 0
beta_R = 500
beta_D = 500
R_crit = 10000
D_crit = 22000
epsilon_U = 0.01
epsilon_Q = 0.01
mu_sa_max = 0.7
mu_sp_max = 0.7
U_crit = 2000
Q_crit = 500
Q_shift = -1500
U_shift = -1500
N_half = 1000
S_half = 2000
alpha_2 = 100

omega_R = log(((R_max - R_min)/epsilon_R) - 1)/R_crit
omega_D = log(((D_max - D_min)/epsilon_D) - 1)/D_crit

omega_U = log((mu_sa_max/0.01) - 1)/U_crit
omega_Q = log((mu_sp_max/0.01) - 1)/Q_crit


#=
x[1] = H
x[2] = N
x[3] = P
x[4] = A
x[5] = S
x[6] = Q
x[7] = U
=#
@taylorize function biomodel7d!(dx, x, p, t)
    I = x[3] + S_n * x[2] - S_a * x[4]
    Rt = R_max - (R_max - R_min)/(1 + exp.(omega_R * I))
    Dt = D_max - (D_max - D_min)/(1 + exp.(omega_D * I))
    omega_U = log((mu_sa_max/epsilon_U) - 1)/U_crit * I
    omega_Q = log((mu_sp_max/epsilon_Q) - 1)/Q_crit * I
    mu_sat = mu_sa_max / (1 + exp(omega_U * (I * U_shift)))
    mu_spt = mu_sp_max / (1 + exp(omega_Q * (I * Q_shift)))
    

    dx[1] = (Rt - Dt) * x[1]
    dx[2] = g_n * x[2] * (1 - x[2]/N_inf) - (k_nq * x[6] + k_ns * x[5]) * x[2]/(N_half + x[2])
    dx[3] = S_PS * x[5] + S_PQ * x[6] + S_PH * x[1] - d_p * x[3]
    dx[4] = S_AU * x[7] + S_AS * x[5] + S_AH * x[1] - d_a * x[4]
    dx[5] = Dt * x[1] - ((1000 * (x[2]/(N_half + x[2]))) * (x[5]/(x[5] + S_half))) - d_s * x[5] - (mu_sat + mu_spt) * x[5]
    dx[6] = mu_spt * x[5] - d_q * x[6]
    dx[7] = mu_sat * x[5] - d_u * x[7]
    return dx
end

#=
@taylorize function biomodel7d!(dx, x, p, t)
    P_S_n = x[3] + S_n
    P_S_n_N = P_S_n * x[2]
    S_a_A = S_a * x[4]
    I = P_S_n_N - S_a_A
    R_range = R_max - R_min
    exponent_R = exp.(omega_R * I)
    exponent_R_1 = (1 + exponent_R)
    R_range_exponent = R_range/exponent_R_1
    Rt = R_max - R_range_exponent
    D_range = D_max - D_min
    exponent_D = exp.(omega_D * I)
    exponent_D_1 = (1 + exponent_D)
    D_range_exponent = D_range/exponent_D_1
    Dt = D_max - D_range_exponent
    mu_epsilon_U = mu_sa_max/epsilon_U
    mu_epsilon_U_1 = mu_epsilon_U - 1
    U_crit_I = U_crit * I
    log_mu_epsilon_U_1 = log(mu_epsilon_U_1)
    omega_U = log_mu_epsilon_U_1/U_crit_I
    mu_epsilon_Q = mu_sp_max/epsilon_Q
    mu_epsilon_Q_1 = mu_epsilon_Q - 1
    Q_crit_I = Q_crit * I
    log_mu_epsilon_Q_1 = log(mu_epsilon_Q_1)
    omega_Q = log_mu_epsilon_Q_1/Q_crit_I
    I_U_shift = I * U_shift
    omega_U_I_U_shift = omega_U * I_U_shift
    exp_omega_U_I_U_shift = exp(omega_U_I_U_shift)
    exp_omega_U_I_U_shift_1 = exp_omega_U_I_U_shift + 1
    mu_sat = mu_sa_max / exp_omega_U_I_U_shift_1
    I_Q_shift = I * Q_shift
    omega_Q_I_Q_shift = omega_Q * I_Q_shift
    exp_omega_Q_I_Q_shift = exp(omega_Q_I_Q_shift)
    exp_omega_Q_I_Q_shift_1 = exp_omega_Q_I_Q_shift + 1
    mu_spt = mu_sp_max / exp_omega_Q_I_Q_shift_1
    Rt_Dt = Rt - Dt
    g_n_N = g_n * x[2]
    N_N_inf = x[2]/N_inf
    N_N_inf_1 = 1 - N_N_inf
    g_n_N_N_N_inf_1 = g_n_N * N_N_inf_1
    k_nq_Q = k_nq * x[6]
    k_ns_S = k_ns * x[5]
    k_nq_Q_k_ns_S = k_nq_Q + k_ns_S
    N_Q_S = g_n_N_N_N_inf_1 - k_nq_Q_k_ns_S
    N_decay = N_half + x[2]
    N_decay_2 = x[2] / (N_decay)
    S_PS_S = S_PS * x[5]
    S_PQ_Q = S_PQ * x[6]
    S_PH_H = S_PH * x[1]
    d_p_P = d_p * x[3]
    S_PS_S_Q = S_PS_S + S_PQ_Q
    S_PH_H_d_p = S_PH_H - d_p_P
    
    S_AS_S = S_AS * x[5]
    S_AU_U = S_AU * x[4]
    S_AH_H = S_AH * x[1]
    d_a_A = d_a * x[4]
    S_AS_S_U = S_AS_S + S_AU_U
    S_AH_H_d_A = S_AH_H - d_a_A

    mu_spt_S = mu_spt * x[5]
    d_Q = d_q * x[6]

    mu_sat_S = mu_sat * x[5]
    d_U = d_u * x[7]

    mu = mu_sat + mu_spt
    mu_S = mu * x[5]
    d_S = d_s * x[5]
    d_S_mu_S = d_S - mu_S
    S_decay = x[5] + S_half
    S_decay_2 = x[5]/ S_decay
    N_S_decay = N_decay_2 * S_decay_2
    N_S_decay_1000 = 1000 * N_S_decay
    N_S_decay_1000_d_S_mu_S = N_S_decay_1000 - d_S_mu_S
    D_T = Dt * x[1]


    dx[1] = Rt_Dt * x[1]
    dx[2] = N_Q_S * N_decay_2
    dx[3] = S_PS_S_Q + S_PH_H_d_p
    dx[4] = S_AS_S_U + S_AH_H_d_A
    dx[5] = D_T - N_S_decay_1000_d_S_mu_S
    dx[6] = mu_spt_S - d_Q
    dx[7] = mu_sat_S - d_U
    return dx
end=#

#X0 = Hyperrectangle(low = [500.0, 0.0, 300, 1000, 5000.0, 0.0, 0.0], high = [1000, 1000, 600, 1160, 6600.0, 200, 200])
#X0 = Hyperrectangle(low = fill(0.9, 7), high = fill(1.5, 7))
X0 = Hyperrectangle(low = [2000.0, 500.0, 100.0, 100.0, 100.0, 100.0, 100.0], high = [2000.1, 500.1, 100.1, 100.1, 100.1, 100.1, 100.1])
function model(X0)
    S = @system(x' = biomodel7d!(x), dim:7)
    return IVP(S, X0)
end


prob = model(X0)
sol = solve(prob, T=100.0);
#sol = solve(prob, AutoVern7(Rodas5()), T=100.0)


H = plot(sol, vars=(0,1), xlabel="time", ylabel="HSPCs")
N = plot(sol, vars=(0,2), xlabel="time", ylabel="Pathogens")
P = plot(sol, vars=(0,3), xlabel="time", ylabel="Pro-Inflammatory")
A = plot(sol, vars=(0,4), xlabel="time", ylabel="Anti-Inflammatory")
S = plot(sol, vars=(0,5), xlabel="time", ylabel="Stable")
Q = plot(sol, vars=(0,6), xlabel="time", ylabel="Active")
U = plot(sol, vars=(0,7), xlabel="time", ylabel="Inactive")
plot(H, N, A, P, S, Q, U, layout = (3, 3), w = 1)

#plot(sol, vars=(0,1), xlabel="time", ylabel="HSPCs")
