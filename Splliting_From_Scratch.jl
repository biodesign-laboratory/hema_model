using IntervalArithmetic
using LazySets
using LinearAlgebra
using Plots
#=
This is an algorithm to run reachability analysis on a hybrid system. The general premise is that you have a list of queues,
each a vector (X, m, t) where X is the set of states at time t and m is the mode the set is being propogated under. X is 
currently just a LazySets.Interval data type. When X intersects with the guard, the part that intersects with the guard goes
into a new queue at the time t that it intersects and it switches to the other mode m. The part that doesn't intersect
is also put into a queue at that time t, but with the mode m it had been running at. These two queues are then added to the
list of queues (i.e. the major_queue). The process then restarts, running each of these queues one after the other, until 
the time T has been reached for all of them.

Problems:
The set propogation is done via linear maps and translations of the intervals. Unsure if this is the most accurate way to do this
The major_queue is full of of queues all at different times t. May be more efficient to attempt to combine these queues
when the time t is the same.
Right now, the script is producing a set of intervals where the upperbound is simply a different magnitude of the time.
First queue is coming back undefined
=#

function plot_res(res)
    #initialize the plot, set of times, set of upper_bounds, and set of lower_bounds
    p = plot()
    times = Float64[]
    upper_bound = Float64[]
    lower_bound = Float64[]

    #fill the sets with each value from the result vector
    for i in eachindex(res)
        interval, time = res[i][1]
        push!(times, time)
        push!(lower_bound, min(interval))
        push!(upper_bound, max(interval))
    end

    #plot the upper and lower sets vs time and connect the two bounds
    plot!(p, lower_bound, fillrange = upper_bound, c = 1)
    return p
end

function run_reach(δ, queue, T, guard)
    #initialize the two queus that will compose of the intervals below and equal to or above the initial interval
    queue_1 = Tuple{LazySets.Interval{Float64}, Int64, Float64}[]
    queue_2 = Tuple{LazySets.Interval{Float64}, Int64, Float64}[]
    res = Tuple{LazySet, Float64}[]
    init, loc, t = pop!(queue)

    #run the initial interval in the continuous function
    R = reach_continuous(loc, init, δ, T-t)
    for i in eachindex(R)
        #Takes each set from t -> T in the continuous state and adds it to S
        #Checks each set to see if it intersects the guard
        S, t = R[i]
        
        #Push the specific set of states S at time t to the result
        push!(res, (S, t))
        if !isdisjoint(S, guard)
            new_t = t + δ * (i + 1)
           
            #intersection of S and guard and splits it then concatenates it
            C = intersection(S, guard)
            C_lower_bound = LazySets.Interval(min(C),10^10)
            C_upper_bound = LazySets.Interval(-10^10,max(C))
            S_lower_bound = S - C_lower_bound
            S_middle = C
            S_upper_bound = S - C_upper_bound
            U = LazySets.Interval(min(S_middle), max(S_upper_bound))
            L = S_lower_bound
            push!(queue_1, (L, 1, new_t))
            push!(queue_2, (U, 2, new_t))
            break
        end
    end
    return (queue_1, queue_2, res)
end

function reach_continuous(loc, init, δ, T)
    #δ is the time step
    #loc in the mode or "location"
    #T is the max time
    #init is is the input set of states
    #Returns a set of states at time t for each T/δ until T is reached

    #discretize the system
    N = floor(Int, T/δ)
    D = [-0.01]
    r = [1.0]

    #Makes a list of H_stable values from the min to the max with 0.01 intervals
    H_stable_temp = Float64[]
    H_stable = init

    #Set the mode
    if loc == 1
        E = 3;
    elseif loc == 2
        E = 2;
    end
    
    # preallocate array
    # R = Vector{Integer}(undef, N)
    R = Vector{Tuple{LazySet, Float64}}(undef, N)
    if N == 0
        return R
    end

    #R[1] = Vector{Tuple{LazySets.MinkowskiSum, Float64}}
    R[1] = (init, 0.0)
    #R[1] = init
    for i in 2:N
        #Translate and scale our initial set based on our mode and system dynamics
        translation_vector = r + D
        dHdt = LazySets.translate(H_stable, translation_vector)
        dHdt = dHdt * E
        
        #Take the derivatives and do a linear approximation to get the set of values at that time
        dHdt = dHdt * N
        translation_vector = translation_vector * i * E
        time = i * δ
        H_stable = LazySets.translate(H_stable, translation_vector)
        R[i] = (H_stable, time)
    end
    return R
end






#Main
#run the major_queue over and over again until the time limit is reached
#declare the major queue and the queue and the res
#queue is a vector (array) of tuples == (Interval, Integer, Float) for (initial interval, loc, time)
#major_queue is an array of vectors of tuples of the same type listed above
queue = Vector{Tuple{LazySets.Interval, Integer, Float64}}(undef, 1)
major_queue = Vector{Tuple{LazySets.Interval, Integer, Float64}}[]
#major_queue = Array{Tuple{LazySets.Interval, Integer, Float64}}(undef, 2, 2)
res = Vector{Tuple{LazySet, Float64}}[]

#Time step, overall time, guard, starting interval
δ = 0.01
T = 4.
guard = LazySets.Interval(0.5, 1.0);
init = LazySets.Interval(0.0, 5.0);

#initialize first queue with the initial interval, mode, and time
queue = (init, 1, 0.0);
push!(major_queue, [queue])

while !isempty(major_queue)
    #Takes latest entry in the major_queue
    queue = pop!(major_queue)
    init, loc, t = queue[1]

    #only run if the time has not exceeded the max time
    if t < 3.9
        queue_immediate = run_reach(δ, queue, T, guard)
    else
        break
    end
    
    #Add the resulting set of states and their times to the result tuple
    push!(res, queue_immediate[3])
    
    #add the two new queues to the major queue
    push!(major_queue, queue_immediate[1])
    push!(major_queue, queue_immediate[2])
end


plot_res(res)
