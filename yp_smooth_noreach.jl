#Smooth model, pure ODE, no reachability.

using Plots, DifferentialEquations

# parameters
K_h = 3
d_H = 0.05
theta_N = 2000
theta_K = 5000
tauQ = 0.5
tauU = 0.5

d_SCSF = 0.3
d_S = 0.8
d_Q = 0.9
d_U = 0.8
d_P = 0.95
d_A = 0.95

g_n = 0.1
N_inf = 20000000
N_half = 500
S_PH = 3
S_PS = 1
S_PQ = 5
S_AU = 7
S_AH = 3
S_AS = 1
s_scsf = 10000
S_KD = 1

k_sn = 3
k_nq = 10
k_ns = 0.5
R_KU = 100
I_crit = 0.2
K_crit = 10000
k = 1
A_crit = 3

H_0 = 1000
N_0 = 0.0
P_0 = 600
A_0 = 1160

# exponent
kk = 1


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
    Hq=x[1];Hm=x[2];N=x[3];S=x[4];Q=x[5];U=x[6];SCSF=x[7];P=x[8];A=x[9];K=x[10]
    
    I = (P * ((N^kk/(N^kk + theta_N^kk)) + 0.25) * ((0.5 * (K/(theta_K + K))) + 1)) / ((A * ((theta_N^kk/(N^kk + theta_N^kk)) + 0.25) * ((0.75 * (K/((theta_K) + K))) + 1)) + (P * ((N^kk/(N^kk + theta_N^kk)) + 0.25) * ((0.5 * (K/((theta_K) + K))) + 1)))
       #(P * ((N^kk)/(theta_N^kk + N^kk) + 0.25) * ((0.5 * (K)/(theta_K + K) + 1))) / ((A * (theta_N^kk/(theta_N^kk + N^kk) + 0.25) * (0.75  *  (K)/(theta_K + K) + 1)) + (P * ((N^kk)/(theta_N^kk + N^kk) + 0.25) * ((0.5  *  (K)/(theta_K + K) + 1))))

    println(I)
    
    nu_Q = 0.2 * ((Hq * P)/(0.2 * Hq + P)) * I
    nu_M = 0.2 * ((Hm * A)/(0.2 * Hm + A)) * (1 - I)
    
    beta = I / ((I_crit + 10) + I)
    
    DI = (Hm*P)/(Hm + 5*P) * I
    DQB = ((1/3) * tauQ * I)/(I_crit + I)
    DUB = ((1/3) * tauU ) * ((1/I) / (A_crit + (1/I)))

    dHq = (nu_M - nu_Q) + (2 * Hq * (1 - (K_h * Hq)/linearlike(SCSF, 0.1, 0.00001)))
    dHm = nu_Q - DI*beta - d_H*Hm - nu_M
    dN = g_n * Hm * (1 - Hm/N_inf) - (k_nq * Q + k_ns * S) * (Hm/(N_half + N))
    dS = (DI * (1 - beta) + 2 * DI * beta) - DQB * S - DUB * S - d_S * S - (k_sn * N * (S/(k_sn * N + S)))
    dQ = DQB * S - d_Q * Q * (1 - 0.5 * I / (2 + I) + 0.5)
    dU = DUB * S - d_U * U
    dSCSF = s_scsf * (K_crit /(K_crit + K)) - d_SCSF * SCSF
    dP = (S_PH * Hm + S_PQ * Q) * (0.5 * I / (2+ I) + 0.5) + S_PS * S - d_P * P
    dA = (S_AH * Hm + S_AU * U) * (0.5 * I / (2+ I) + 0.5) + S_AS * S - d_A * A
    dK = S_KD * (k_sn * N * (S/(k_sn * N + S))) - (R_KU * U * (K/K_crit * K))

    dx[1]=dHq;dx[2]=dHm;dx[3]=dN;dx[4]=dS;dx[5]=dQ;dx[6]=dU;dx[7]=dSCSF;dx[8]=dP;dx[9]=dA;dx[10]=dK
    return dx
end

dx_test = Any[0,0,0,0,0,0,0,0,0,0]
x_test = Any[1,2,3,4,5,6,7,8,9,10]

# compare output to Bryan's RHS
println(biomodel(dx_test,x_test,0,0))

