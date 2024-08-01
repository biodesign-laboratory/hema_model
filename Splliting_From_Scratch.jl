using IntervalArithmetic
using LazySets
using LinearAlgebra
using Plots

#= This version has Pathogens as the interval,
everything else is scalar=#

function plot_res(res)
    #initialize the plot, set of times, set of upper_bounds, and set of lower_bounds
    p = plot()
    times = Float64[]
    upper_bound = Float64[]
    lower_bound = Float64[]

    #fill the sets with each value from the result vector
    for i in eachindex(res)
        interval, time = res[i]
        push!(times, time)
        push!(lower_bound, min(interval))
        push!(upper_bound, max(interval))
        #plot!(p, [time, time], [min(interval), max(interval)], fillrange=lower_bound, fillalpha=0.3, label="", color="blue", xlims=[0,0.5])
        plot!(p, [time, time], [min(interval), max(interval)], fillrange=lower_bound, fillalpha=0.3, label="", color="blue")
    end


    #plot the upper and lower sets vs time and connect the two bounds
    #plot the intervals vs time and connect the bounds
    #plot!(p, times, lower_bound, ribbon=(upper_bound), fillalpha=0.3, label="Interval bounds")
    xlabel!("Time")
    ylabel!("Interval")
    title!("Intervals vs Time")
    #print(res)
    return p
end

function run_reach(δ, local_queue, T, guard)
    #initialize the two queus that will compose of the intervals below and equal to or above the initial interval
    global rates
    NaN_bool = false
    res = Vector{Tuple{LazySet, Float64}}(undef, 1)
    queue_1 = (LazySets.Interval(0.0, 1.0), 1, 0.0)
    queue_2 = (LazySets.Interval(0.0, 1.0), 1, 0.0)
    init, loc, t = last(local_queue)

    #run the initial interval in the continuous function
    R, rate_times = reach_continuous(loc, init, δ, T, t)
    res[1] = (init, t)
    for i in eachindex(R)
        #Takes each set from t -> T in the continuous state and adds it to S
        #Checks each set to see if it intersects the guard
        S, t = R[i]

        if low(S)[1] < 0 && high(S)[1] > 0
            S = LazySets.Interval(0.0, high(U)[1])
        elseif low(S)[1] < 0 && high(S)[1] < 0
            S = LazySets.Interval(0.0, 0.0)
        end
        
        #Push the specific set of states S at time t to the result
        push!(res, (S, t))

        #making a new S_temp to reflect that we aren't checking for P, but for a translated version of P
        S_temp = S
        try
            scaling_factor = get(rates, :S_n, undef) + rate_times[2][i]
            println("scaling factor: ", scaling_factor)
            S_temp = LazySets.Interval(scaling_factor * min(S), scaling_factor * max(S))
            translation_vector = [get(rates, :S_a, undef) * rate_times[3][i]]
            #println("translation: ", translation_vector)
            S_temp = LazySets.translate(S_temp, translation_vector)
            println(S_temp)
        catch e
            println(e)
            break
        end
        if !isdisjoint(S_temp, guard)
            println("split")
            new_t = R[i][2] + δ
            L = LazySets.Interval(0.0, 1.0)
            U = LazySets.Interval(0.0, 1.0)

            if low(S_temp) > low(guard)
                U = S_temp
          
                U = scale(1/scaling_factor, U)
                U = LazySets.translate(U, -translation_vector)

                if low(U)[1] < 0 && high(U)[1] > 0
                    U = LazySets.Interval(0.0, high(U)[1])
                elseif low(U)[1] < 0 && high(U)[1] < 0
                    U = LazySets.Interval(0.0, 0.0)
                end

                queue_2 = (U, 2, new_t)
                queue_1 = queue_2
            elseif low(S_temp) < low(guard) && high(S_temp) != high(guard)
                L = LazySets.Interval(low(S_temp)[1], prevfloat(low(guard)[1]))
                U = LazySets.Interval(prevfloat(low(guard)[1]), nextfloat(high(S_temp)[1]))

                U = scale(1/scaling_factor, U)
                U = LazySets.translate(U, -translation_vector)
                L = scale(1/scaling_factor, L)
                L = LazySets.translate(L, -translation_vector)

                if low(U)[1] < 0 && high(U)[1] > 0
                    U = LazySets.Interval(0.0, high(U)[1])
                elseif low(U)[1] < 0 && high(U)[1] < 0
                    U = LazySets.Interval(0.0, 0.0)
                end
                if low(L)[1] < 0 && high(L)[1] > 0
                    L = LazySets.Interval(0.0, high(L)[1])
                elseif low(L)[1] < 0 && high(L)[1] < 0
                    L = LazySets.Interval(0.0, 0.0)
                end

                queue_2 = (U, 2, new_t)
                queue_1 = (L, 1, new_t)
            end

            #Return rate values to where they were at the time S hit the guard
            rates[:H] = rate_times[1][i]
            rates[:P] = rate_times[2][i]
            rates[:A] = rate_times[3][i]
            nodes[:S] = rate_times[4][i]
            nodes[:Q] = rate_times[5][i]
            nodes[:U] = rate_times[6][i]
            NaN_bool = check_NaN(rates)
            break
        end
        if i == length(R)
            break
        end
    end
    return (queue_1, queue_2, res, NaN_bool)
end

function reach_continuous(loc, init, δ, T_max, t_0)
    #δ is the time step
    #loc is the mode or "location"
    #T is the max time
    #init is is the input set of states
    #Returns a set of states at time t for each T/δ until T is reached
    global rates
    global nodes

    #discretize the system
    τ = T_max - t_0
    N = floor(Int, τ/δ)

    #Makes a list of P values from the min to the max with 0.01 intervals
    N_interval = init

    #initialize time_rate bookkeeping vectors
    rate_times = []

    H_times = []
    P_times = []
    A_times = []
    S_times = []
    Q_times = []
    U_times = []

    push!(H_times, get(rates, :H, undef))
    push!(P_times, get(rates, :P, undef))
    push!(A_times, get(rates, :A, undef))
    push!(S_times, get(nodes, :S, undef))
    push!(Q_times, get(nodes, :Q, undef))
    push!(U_times, get(nodes, :U, undef))

    #Set the mode
    if loc == 1
        E = 1
    elseif loc == 2
        E = 2 - (2/(1+exp(-1 * get(rates, :y, undef) * get(rates, :P, undef) + get(rates, :S_n, undef) * (LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef) - get(rates, :P_crit, undef)))))
    end

    if E < 0
        E = 0
    end
    
    #preallocate array
    #if N == 0 is the stopping criteria
    R = Vector{Tuple{LazySet, Float64}}(undef, 1)
    R[1] = (init, t_0)
    if N < δ
        return R
    end

    for i in 2:N
        #Time at each i
        time = (i * δ) + t_0

        #Translate and scale our initial set based on our mode and system dynamics
        #add up the translation vector i times and multiply by E similar to the dynamics above
        Renewal = 0
        if (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= 0
            Renewal = 0
        elseif (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= get(rates, :P_crit, undef)
            Renewal = (0.35 * ((get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)/get(rates, :P_crit, undef)))) * get(rates, :H, undef)
        elseif (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) > get(rates, :P_crit, undef)
            Renewal = (0.35 * get(rates, :H, undef))
        end

        D = 0
        if (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= 0
            D = 0
        elseif (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= get(rates, :P_crit, undef)
            D = (0.35 * ((get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef))/get(rates, :P_crit, undef))) * get(rates, :H, undef)
        elseif (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) > get(rates, :P_crit, undef)
            D = (0.35 * get(rates, :H, undef))
        end

        dE_star = 0
        if (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= (1.5 * get(rates, :P_crit, undef)) || get(rates, :E_star, undef) <= 0.1
            dE_star = 0
        else
            dE_star = -0.0005 * (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef) - 1.5 * get(rates, :P_crit, undef))
        end

        rates[:E_star] = rates[:E_star] + (dE_star * time)

        if get(rates, :E_star, undef) * get(rates, :H_0, undef) > get(rates, :H, undef)
            H_stable = 0.05 * get(rates, :H, undef) * (1 - (get(rates, :H, undef)/(get(rates, :E_star, undef) * get(rates, :H_0, undef))))
        else
            H_stable = 0
        end

        A_l = (get(rates, :k_tn, undef) * LazySets.center(N_interval)[1]) * ((-2 / (1 + exp(get(rates, :w, undef) * get(nodes, :S, undef)))) + 1)

        mu_SP = 0
        if (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= get(rates, :Active_start, undef)
            mu_SP = 0
        elseif (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= get(rates, :Active_crit, undef)
            mu_SP = 0.5 * ((get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef) - get(rates, :Active_start, undef)) / get(rates, :Active_crit, undef)) * get(nodes, :S, undef)
        else
            mu_SP = 0.5 * get(nodes, :S, undef)
        end

        mu_SA = 0
        if (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= get(rates, :Immune_start, undef)
            mu_SA = 0
        elseif (get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef)) <= get(rates, :Immune_crit, undef)
            mu_SA = 0.2 * ((get(rates, :P, undef) + get(rates, :S_n, undef) * LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef) - get(rates, :Immune_start, undef)) / get(rates, :Immune_crit, undef)) * get(nodes, :S, undef)
        else
            mu_SA = 0.2 * get(nodes, :S, undef)
        end        

        D_P = get(rates, :d_p, undef) * get(rates, :P, undef)

        if get(rates, :A, undef) > 0
            D_A = get(rates, :d_a, undef) * get(rates, :A, undef)
        else
            D_A = 0
        end

        D_S = get(rates, :d_s, undef) * get(nodes, :S, undef)
        D_Q = get(rates, :d_q, undef) * get(nodes, :Q, undef)
        D_U = get(rates, :d_u, undef) * get(nodes, :U, undef)

        #Node rate calculations
        dHdt = E * (Renewal + H_stable) - D
        dNdt = get(rates, :g_N, undef) * LazySets.center(N_interval)[1] - ((get(rates, :k_nq, undef) * get(nodes, :Q, undef) + (get(rates, :k_ns, undef) * get(nodes, :S, undef))) * (1 - ((LazySets.center(N_interval)[1]/get(rates, :N_inf, undef)))))
        dPdt = (get(rates, :S_PS, undef) * get(nodes, :S, undef)) + (get(rates, :S_PQ, undef) * get(nodes, :Q, undef)) + (get(rates, :S_PH, undef) * get(rates, :H, undef)) - D_P
        dAdt = (get(rates, :S_AU, undef) * get(nodes, :U, undef)) + (get(rates, :S_AS, undef) * get(nodes, :S, rates)) + (get(rates, :S_AH, undef) * get(rates, :H, undef)) - D_A
        dSdt = D - A_l - D_S - mu_SA - mu_SP
        dQdt = mu_SP - D_Q
        dUdt = mu_SA - D_U

        #Linear approximation
        time = (i * δ) + t_0
        translation_vector = [dNdt] * time

        rates[:H] = (dHdt * time) + get(rates, :H, undef)
        rates[:P] = (dPdt * time) + get(rates, :P, undef)
        rates[:A] = (dAdt * time) + get(rates, :A, undef)
        nodes[:S] = (dSdt * time) + get(nodes, :S, undef)
        nodes[:Q] = (dQdt * time) + get(nodes, :Q, undef)
        nodes[:U] = (dUdt * time) + get(nodes, :U, undef)

        
        #Initialize these and append
        push!(H_times, get(rates, :H, undef))
        push!(P_times, get(rates, :P, undef))
        push!(A_times, get(rates, :A, undef))
        push!(S_times, get(nodes, :S, undef))
        push!(Q_times, get(nodes, :Q, undef))
        push!(U_times, get(nodes, :U, undef))
        
        rate_times = [H_times, P_times, A_times, S_times, Q_times, U_times]

        try
            N_interval = LazySets.translate(N_interval, translation_vector)
        catch err
            if isa(err, AssertionError)
                println("AssertionError: Attempted to access an index that is out of bounds. Error details: ", err)
                push!(R, (N_interval, time))
                break
            else
                rethrow(err)  # Rethrow the error if it is not a BoundsError
                break
            end
        end

        #Hybrid dynamic
        if loc == 1
            E = 1
        elseif loc == 2
            E = 2 - (2/(1+exp(-1 * get(rates, :y, undef) * get(rates, :P, undef) + get(rates, :S_n, undef) * (LazySets.center(N_interval)[1] - get(rates, :S_a, undef) * get(rates, :A, undef) - get(rates, :P_crit, undef)))))
        end
    
        if E < 0
            E = 0
        end

        if high(N_interval)[1] > 1.0E+10
            break
        elseif low(N_interval)[1] == 0.0
            break
        end
        if low(N_interval)[1] < 0 && high(N_interval)[1] > 0
            N_interval = LazySets.Interval(0.0, high(N_interval)[1])
        elseif low(N_interval)[1] < 0 && high(N_interval)[1] < 0
            N_interval = LazySets.Interval(0.0, 0.0)
        end
        push!(R, (N_interval, time))
        print(E)
    end
    return R, rate_times
end

function check_NaN(dict)
    NaN_bool = false
    for (key, value) in dict
        if isnan(value)
            NaN_bool = true
        end
    end
    return NaN_bool
end

function main(queue, res, guard, T, δ)
    while !isempty(queue)
        #Takes last entry in the queue
        init, loc, t = last(queue)
    
        #only run if the time has not exceeded the max time
        if t < 3.9
            queue_immediate = run_reach(δ, queue, T, guard)
        else
            break
        end
        
        pop!(queue)
    
        #Add the resulting set of states and their times to the result tuple
        for i in eachindex(queue_immediate[3])
            push!(res, queue_immediate[3][i])
        end
        
        #add the two new queues to the major queue
        #checks if the queues are equal from the logic in the function
        if queue_immediate[1][1] != queue_immediate[2][1]
            push!(queue, queue_immediate[1])
            push!(queue, queue_immediate[2])
        elseif queue_immediate[1][1] != LazySets.Interval(0.0, 1.0) && queue_immediate[2][1] != LazySets.Interval(0.0, 1.0)
            push!(queue, queue_immediate[2])
        end

        if queue_immediate[4] == true
            break
        end
    end
end


#Set rates
rates = Dict(
        :w => 0.0005,
        :P_crit => 7000,
        :A_crit => 1000,
        :t_half_leukocytes => 7,
        :t_half_P => 4.1375,
        :t_half_A => 7,
        :t_double_P => 15,
        :gamma => 0.0000001,
        :S_a => 1,
        :S_n => 1,
        :N_inf => 20000000,
        :S_PQ => 0.33,
        :S_PH => 0.01,
        :S_PS => 0.02,
        :S_AS => 0.04,
        :S_AH => 0.01,
        :S_AU => 0.33,
        :theta_ps => 10000000,
        :theta_ar => 10000000,
        :theta_AS => 10000000,
        :theta_UP => 10000000,
        :Immune_start => 500,
        :Active_start => 500,
        :Immune_crit => 2500,
        :Active_crit => 2500,
        :y => 2500,
        :d_s => 1/70,
        :d_p => 1/4,
        :d_a => 1/4,
        :d_q => 1/4,
        :d_u => 0.05,
        :g_N => 0.2,
        :K_PS => 0.7,
        :K_AQ => 0.7,
        :K_PU => 0.7,
        :K_AS => 0.7,
        :k_nq => 0.85,
        :k_ns => 0.2,
        :k_tn => 0.33,
        :H => 1000,
        :N_func => 0.0,
        :P => 600,
        :A => 1160,
        :T => 7000,
        :a => 0.00,
        :b => 1,
        :e => 0.00,
        :E_star => 1,
        :H_0 => 1000,
        :N_0 => 0.0,
        :P_0 => 600,
        :A_0 => 1160,
        :T_0 => 7000,
        :a_0 => 0.00,
        :b_0 => 1,
        :e_0 => 0.00,
        :E_star_0 => 1
    )
nodes = Dict(
    :Q => get(rates, :a, undef) * get(rates, :T, undef),
    :S => get(rates, :b, undef) * get(rates, :T, undef),
    :U => get(rates, :e, undef) * get(rates, :T, undef)
)

#Time step, overall time, guard, starting interval, queue and results
queue = Vector{Tuple{LazySets.Interval, Integer, Float64}}(undef, 1)
res = Vector{Tuple{LazySet, Float64}}(undef, 1)
δ = 0.01
T = 4.
guard = LazySets.Interval((get(rates, :P_crit, undef) - 5), (get(rates, :P_crit, undef) + 5))
init = LazySets.Interval(500, 1000)

#initialize first queue with the initial interval, mode, and time
queue[1] = (init, 1, 0.0)
init, loc, t = queue[1]
res[1] = (init, t)

main(queue, res, guard, T, δ)
plot_res(res)
