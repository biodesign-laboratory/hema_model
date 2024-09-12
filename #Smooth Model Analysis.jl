#Smooth Model Analysis
using ReachabilityAnalysis, ModelingToolkit, Plots

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
    S_n = 1
    N_inf = 20000000
    S_PQ = 0.33
    S_PH = 0.01
    S_PS = 0.02
    S_AS = 0.04
    S_AH = 0.01
    S_AU = 0.33
    theta_ps = 10000000
    theta_ar = 10000000
    theta_AS = 10000000
    theta_UP = 10000000
    Immune_start = 500
    Active_start = 500
    Immune_crit = 2500
    Active_crit = 2500
    y = 2500
    d_s = 1/70
    d_p = 1/4
    d_a = 1/4
    d_q = 1/4
    d_u = 0.05
    g_n = 0.2
    K_PS = 0.7
    K_AQ = 0.7
    K_PU = 0.7
    K_AS = 0.7
    k_nq = 0.85
    k_ns = 0.2
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

#Variables for H
epsilon_R = 0.001
epsilon_D = 0.001
R_max = 5000
R_min = 0
D_max = 10500
D_min = 0
beta_R = 500
beta_D = 500
R_crit = 200
D_crit = 200

#Variables for N
#N_half = (k_nq * x[6] + k_ns * x[5])/2

#Variables for S Q U
#S_half = (k_nq * x[6] + k_ns * x[5])/2
#A_l = k_tn * N * (S / (S_half + S))
epsilon_U = 0.01
epsilon_Q = 0.01
mu_sa_max = 10000
mu_sp_max = 10000
U_crit = 1000
Q_crit = 1000
#omega_U = log((mu_sa_max/epsilon_U) - 1)/U_crit * I
#omega_Q = log((mu_sp_max/epsilon_Q) - 1)/Q_crit * I
#mu_sa = mu_sa_max * 1/(1 + exp(omega_U))
#mu_sp = mu_sp_max * 1/(1 + exp(omega_Q))

#Prep work for H
#=I = P + S_n * N - S_a * A
omega_R = log(((R_max - R_min)/epsilon_R) - 1)/R_crit * (I - beta_R)
omega_D = log(((D_max - D_min)/epsilon_D) - 1)/D_crit * (I - beta_D)
R = R_max - (R_max - R_min)/(1 + exp.(omega_R))
D = D_max - (D_max - D_min)/(1 + exp.(omega_D))
=#

#Prep work for N
#x[2] = g_n * x[2](1 - N_inf) - (k_nq * x[6] + k_ns * x[5]) * x[2]/(N_half + x[2])


#Prep work for S Q U
#dx[5] = D - A_l - D_s - mu_sa - mu_sp
#dQdt = mu_SP - D_Q
#dUdt = mu_SA - D_U

#Prep work for P and A
#D_p = d_p * x[3]
#dx[3] = S_PS * x[5] + S_PQ * x[6] + S_PH * x[1] - D_p
#D_a = d_a * A
#dx[4] = S_AU * x[7] + S_AS * x[5] + S_AH * x[1] - D_a

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
    omega_R = log(((R_max - R_min)/epsilon_R) - 1)/R_crit * (I - beta_R)
    omega_D = log(((D_max - D_min)/epsilon_D) - 1)/D_crit * (I - beta_D)
    R = R_max - (R_max - R_min)/(1 + exp.(omega_R))
    D = D_max - (D_max - D_min)/(1 + exp.(omega_D))

    N_half = 1000
    S_half = 1000
    A_l = k_tn * x[2] * (x[5] / (S_half + x[5]))
    omega_U = log((mu_sa_max/epsilon_U) - 1)/U_crit * I
    omega_Q = log((mu_sp_max/epsilon_Q) - 1)/Q_crit * I
    mu_sa = mu_sa_max * 1/(1 + exp(omega_U))
    mu_sp = mu_sp_max * 1/(1 + exp(omega_Q))
    D_s = d_s * x[5]
    D_q = d_q * x[6]
    D_u = d_u * x[7]
    D_p = d_p * x[3]
    D_a = d_a * x[4]


    dx[1] = (R - D) * H
    dx[2] = g_n * x[2] * (1 - N_inf) - (k_nq * x[6] + k_ns * x[5]) * x[2]/(N_half + x[2])
    dx[3] = S_PS * x[5] + S_PQ * x[6] + S_PH * x[1] - D_p
    dx[4] = S_AU * x[7] + S_AS * x[5] + S_AH * x[1] - D_a
    dx[5] = D - A_l - D_s - mu_sa - mu_sp
    dx[6] = mu_sp - D_q
    dx[7] = mu_sa - D_u
    return dx
end

#X0 = Hyperrectangle(low = [500, 0.0, 300, 1000, 5000.0, 0.0, 0.0], high = [1000, 1000, 600, 1160, 6600.0, 200, 200])
X0 = Hyperrectangle(low = fill(0.9999, 7), high = fill(1.0001, 7))
function model(X0)
    S = @system(x' = biomodel7d!(x), dim:7)
    return IVP(S, X0)
end


prob = model(X0)
sol = solve(prob, T=1.0, abstol = 1e-5)

plot(sol, vars=(0,1))