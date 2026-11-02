# Find the upper nearest of τ element of time 
"""
    nearest_time(τ::Vector{Float64}, time::Vector{Float64})

τ: A vector of Float64 values representing the times to be compared.
time: A vector of Float64 values representing the reference times.

Returns the sequence of indices of the first elements of time that are greater than or equal to each element of τ.

τ and time need to be sorted in ascending order.
"""
function nearest_time(τ::Vector{Float64}, time::Vector{Float64})

    #creation of variables for the next step
    n = length(τ)
    time_stop_index = 1
    time_stop_indices = Vector{Int}(undef, n)

    #approximate τ by the the first index of time above τ[i]    
    for i ∈ 1:n
        while τ[i] > time[time_stop_index]
            time_stop_index += 1
        end
        time_stop_indices[i] = time_stop_index 
    end

    return time_stop_indices
end

"""
    nearest_time(τ::Float64, time::Vector{Float64})

τ: A vector of Float64 values representing the times to be compared.
time: A vector of Float64 values representing the reference times.

Returns the index of the first element of time that is greater than or equal to τ.

time need to be sorted in ascending order.
"""
nearest_time(τ::Float64, time::Vector{Float64}) = nearest_time([τ], time)[1]

"""
    nearest_time(τ::Float64, time::Vector{Float64})

τ: A vector of Float64 values representing the times to be compared.
time: A vector of Float64 values representing the reference times.

Returns the index of the first element of time that is greater than or equal to τ.

time need to be sorted in ascending order.
"""
nearest_time(time::Vector{Float64}, τ::Float64) = nearest_time([τ], time)[1]