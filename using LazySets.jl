using LazySets

function split_intervals(S::LazySets.Interval, g::LazySets.Interval)
    intervals = LazySets.Interval{Float64}[]

    # Part below g
    if low(S) < low(g)
        push!(intervals, LazySets.Interval(low(S)[1], prevfloat(low(g)[1])))
    end

    # Part above g
    if high(S) > high(g)
        push!(intervals, LazySets.Interval(nextfloat(high(g)[1]), high(S)[1]))
    end

    return intervals
end

# Example usage
S = LazySets.Interval(0.5, 1.0)
g = LazySets.Interval(0.8, 0.9)
#x = LazySets.Interval(low(g)[1], high(g)[1])
#h = x + g
#print(x)
result = split_intervals(S, g)
println(result)  # Output: [0.5, 0.7999999999999999], [0.9, 1.0]
