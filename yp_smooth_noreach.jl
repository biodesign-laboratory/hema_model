#Smooth model, pure ODE, no reachability.

using Plots, DifferentialEquations

# Extract rates into individual variables
    A_crit = 3
    N_inf = 20000000
    N_half = 500
    S_PQ = 5
    S_PH = 3
    S_PS = 1
    S_AS = 1
    S_AH = 3
    S_AU = 7
    
    d_s = 0.8
    d_p = 0.95
    d_a = 0.95
    d_q = 0.9
    d_u = 0.8
    g_n = 0.1
    
    k_nq = 10
    k_ns = 0.5
    k_sn = 3

    H_0 = 1000
    N_0 = 0.0
    P_0 = 600
    A_0 = 1160

#Variables
theta_N = 2000
theta_K = 5000
k = 1
K_h = 3
K_crit = 10000
d_scsf = 0.3
s_scsf = 10000
d_h = 0.05
S_KD = 1
R_KU = 100
kk = 1
I_crit = 0.2
tauQ = 0.5
tauU = 0.5

function linearlike(x,c,e)
    c + (x - c) * (x^2 / (x^2 + e))
end

#=
x[1] = Hq
x[2] = Hm
x[3] = N
x[4] = S
x[5] = Q
x[6] = U
x[7] = SCSF
x[8] = P
x[9] = A
x[10] = K
=#


function biomodel(dx, x, p, t)
    I = (x[8] * ((x[3]^kk/(x[3]^kk + theta_N^kk)) + 0.25) * ((0.5 * (x[10]/((theta_K) + x[10]))) + 1)) / ((x[9] * ((theta_N^kk/(x[3]^kk + theta_N^kk)) + 0.25) * ((0.75 * (x[10]/((theta_K) + x[10]))) + 1)) + (x[8] * ((x[3]^kk/(x[3]^kk + theta_N^kk)) + 0.25) * ((0.5 * (x[10]/((theta_K) + x[10]))) + 1)))
    nu_M = 0.2 * ((x[2] * x[9])/(0.2 * x[2] + x[9])) * (1 - I)
    nu_Q = 0.2 * ((x[1] * x[8])/(0.2 * x[1] + x[8])) * I
    beta = I / ((I_crit + 10) + I)
    D_q = d_q * x[5]
    D_u = d_u * x[6]
    D_h = d_h * x[2]
    D_p = d_p * x[8]
    D_a = d_a * x[9]
    D_s = d_s * x[4]
    DI = (x[2] * (1/5) * x[8])/(x[2] + (1/5) * x[8]) * I
    DQB = ((1/3) * tauQ * I)/(I_crit + I)
    DUB = ((1/3) * tauU ) * ((1/I) / (A_crit + (1/I)))

    dx[1] = (nu_M - nu_Q) + (2 * x[1] * (1 - (K_h * x[1])/linearlike(x[7], 0.1, 0.00001)))
    dx[2] = nu_Q - (DI * beta) - D_h - nu_M
    dx[3] = g_n * x[2] * (1 - x[2]/N_inf) - (k_nq * x[5] + k_ns * x[4]) * (x[2]/(N_half + x[3]))
    dx[4] = (DI * (1 - beta) + 2 * DI * beta) - DQB * x[4] - DUB * x[4] - D_s - (k_sn * x[3] * (x[4]/(k_sn * x[3] + x[4])))
    dx[5] = DQB * x[4] - D_q * (1 - 0.5 * I / (2 + I) + 0.5)
    dx[6] = DUB * x[4] - D_u
    dx[7] = s_scsf * (K_crit /(K_crit + x[10])) - d_scsf * x[7]
    dx[8] = (S_PH * x[2] + S_PQ * x[5]) * (0.5 * I / (2+ I) + 0.5) + S_PS * x[4] - D_p
    dx[9] = (S_AH * x[2] + S_AU * x[6]) * (0.5 * I / (2+ I) + 0.5) + S_AS * x[4] - D_a
    dx[10] = S_KD * (k_sn * x[3] * (x[4]/(k_sn * x[3] + x[4]))) - (R_KU * x[6] * (x[10]/K_crit * x[10]))
    return dx
end

dx_test = Any[1,2,3,4,5,6,7,8,9,10]
x_test = Any[1,2,3,4,5,6,7,8,9,10]

# compare output to Bryan's RHS
println(biomodel(dx_test,x_test,0,0))

