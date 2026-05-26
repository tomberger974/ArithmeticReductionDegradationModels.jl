# The 2 following functions exists to find the upper nearest of τ element of time 
function nearest_time(τ::Float64, time::Vector{Float64})

    #approximate τ by the the first index of time above τ
    time_stop_index = 1
    while τ > time[time_stop_index]
        time_stop_index += 1
    end

    return time_stop_index
end

function nearest_time(τ::Vector{Float64}, time::Vector{Float64})

    #creation of variables for the next step
    n = length(τ)
    time_stop_index = 1
    time_stop_indices = ones(Int, n)

    #approximate τ by the the first index of time above τ[i]    
    for i ∈ 1:n
        while τ[i] > time[time_stop_index]
            time_stop_index += 1
        end
        time_stop_indices[i] = time_stop_index 
    end

    return time_stop_indices
end