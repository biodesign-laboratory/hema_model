#Smooth Model Analysis #2
    using ReachabilityAnalysis, ModelingToolkit, Plots, OrdinaryDiffEq, DifferentialEquations, DiffEqCallbacks

    # Extract rates into individual variables
        w = 0.0005
        P_crit = 7000
        A_crit = 3
        t_half_leukocytes = 7
        t_half_P = 4.1375
        t_half_A = 7
        t_double_P = 15
        gamma = 0.0000001
        S_a = 1
        S_n = 4
        N_inf = 20000000
        N_half = 500
        S_PQ = 5
        S_PH = 3
        S_PS = 1
        S_AS = 1
        S_AH = 3
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
        d_s = 0.8
        d_p = 0.95
        d_a = 0.95
        d_q = 0.9
        d_u = 0.8
        g_n = 0.1
        K_PS = 0.7
        K_AQ = 0.7
        K_PU = 0.7
        K_AS = 0.7
        k_nq = 10
        k_ns = 0.5
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

    @taylorize function biomodel7d!(dx, x, p, t)
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
        DI = (x[2] * (1/5) * x[8])/(x[2] * (1/5) + x[8]) * I
        DQB = ((1/3) * tauQ * I)/(I_crit + I)
        DUB = ((1/3) * tauU  * (1 - I)) / (A_crit + (1 - I))

        dx[1] = (nu_M - nu_Q) + (2 * x[1] * (1 - (K_h * x[1])/linearlike(x[7], 0.1, 0.00001)))
        dx[2] = nu_Q - (DI * beta) - D_h - nu_M
        dx[3] = g_n * x[3] * (1 - x[3]/N_inf) - (k_nq * x[5] + k_ns * x[4]) * (x[3]/(N_half + x[3]))
        dx[4] = (DI * (1 - beta) + 2 * DI * beta) - DQB * x[4] - DUB * x[4] - D_s - (k_sn * x[3] * ((1 - (DQB + DUB + d_s)) * x[4] / (k_sn * x[3] + (1 - (DQB + DUB + d_s)) * x[4])))
        #dx[4] = DI * (1 - beta) + 2 * DI * beta - DQB * x[4] - DUB * x[4] - D_s - (k_sn * x[3] * x[4]/(k_sn * x[3] + x[4]))
        dx[5] = DQB * x[4] - D_q * (1 - 0.5 * I / (2 + I))
        dx[6] = DUB * x[4] - D_u
        dx[7] = s_scsf * (K_crit /(K_crit + x[10])) - d_scsf * x[7]
        dx[8] = (S_PH * x[2] + S_PQ * x[5]) * (0.5 * I / (2 + I) + 0.5) + S_PS * x[4] - D_p
        dx[9] = (S_AH * x[2] + S_AU * x[6]) * (0.5 * I / (2 + I) + 0.5) + S_AS * x[4] - D_a
        dx[10] = S_KD * (k_sn * x[3] * (x[4]/(k_sn * x[3] + x[4]))) - (R_KU * x[6] * (x[10]/K_crit * x[10]))
        return dx
    end

    #X0 = Hyperrectangle(low = fill(0.99, 10), high = fill(1.00, 10))
    X0 = Hyperrectangle(low = [11000.0, 300.0, 0.0, 10.0, 1.0, 0.5, 32000.0, 1000.0, 1000.0, 0.0], high = [11000.1, 300.1, 100.1, 10.1, 1.1, 0.6, 32000.1, 1000.1, 1000.1, 0.1])

    function model(X0)
        S = @system(x' = biomodel7d!(x), dim:10)
        return IVP(S, X0)
    end

   
    prob = model(X0)
    sol = solve(prob, T=50.0);

    try
        H_m = plot(sol, vars=(0,1), xlabel="time", ylabel="H_m")
        H_q = plot(sol, vars=(0,2), xlabel="time", ylabel="H_q")
        N = plot(sol, vars=(0,3), xlabel="time", ylabel="PAMPS")
        S = plot(sol, vars=(0,4), xlabel="time", ylabel="Stable")
        Q = plot(sol, vars=(0,5), xlabel="time", ylabel="Active")
        U = plot(sol, vars=(0,6), xlabel="time", ylabel="Inactive")
        SCSF = plot(sol, vars=(0,7), xlabel="time", ylabel="SCSF")
        P = plot(sol, vars=(0,8), xlabel="time", ylabel="Pro-Inflammatory")
        A = plot(sol, vars=(0,9), xlabel="time", ylabel="Anti-Inflammatory")
        K = plot(sol, vars=(0,10), xlabel="time", ylabel="DAMPS")
        plot(H_m, H_q, N, S, Q, U, SCSF, P, A, K, layout = (4, 3), w = 1)
    catch e
        print(e)
        H_m = plot(sol, vars=(0,1), xlabel="time", ylabel="H_m")
        plot(H_m)
    end
