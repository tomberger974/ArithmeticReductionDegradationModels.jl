#####TRAJECTOIRE D'UN PROCESSUS DE WIENER#####
function nl_wiener_rand(nlwienerard::NLWienerARD, time::Vector{Float64})
    k = length(time)
    μ, σ, α, β = nlwienerard.drift, nlwienerard.dispertion, nlwienerard.power_drift, nlwienerard.power_dispersion

    Δtα = [time[i]^α - time[i-1]^α for i in 2:k]
    Δtβ = [time[i]^β - time[i-1]^β for i in 2:k]
    
    #creation of a sample of the distribution 
    X = rand.(Normal.(μ .* Δtα, σ .* sqrt.(Δtβ)))

    return cumsum(vcat(0., X))
end

# Works only if there is no before nor after maintenance observation
function rand!(nlwienerard::NLWienerARD, degradationdata::DegradationData)

    # Extract model parameters and data
    ρ = nlwienerard.efficiencies
    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    # Create time vector with tags to track element types
    n_deg = nrow(deg)
    n_maint = nrow(maint)
    time = vcat(0., deg.DATE, maint.DATE)
    time_types = vcat(:maint, fill(:deg, n_deg), fill(:maint, n_maint))

    # Sort time and track types
    ordered_time_index = sortperm(time)
    ordered_time = time[ordered_time_index]
    ordered_time_types = time_types[ordered_time_index]

    # Find positions of maintenance and degradation dates in sorted array
    maint_index = findall(==(Symbol(:maint)), ordered_time_types)
    deg_index = findall(==(Symbol(:deg)), ordered_time_types)

    # Generate Wiener process values at all time points
    Y = nl_wiener_rand(nlwienerard, ordered_time)

    # Adjust for maintenance effects
    for i in eachindex(maint_index[1:end-1])
        if typeof(nlwienerard) == NLWienerARD1
                a = Y[maint_index[i + 1]] - Y[maint_index[i]]
        else
                a = Y[maint_index[i + 1]]
        end
        ρi = ρ[maint[i, "TYPE"]]
        Y[maint_index[i+1]:end] .-= ρi*a
    end

    degradationdata.degradations.VALUE .= Y[deg_index]

    return degradationdata
end

# Same function as abosve with just a sample of inspection dates and maintenance dates
function rand(nlwienerard::NLWienerARD, inspection_dates::Vector{Float64}, maint::DataFrame)

    # Create time vector with tags to track element types
    n_deg = length(inspection_dates)
    n_maint = nrow(maint)
    time = vcat(0., inspection_dates, maint.DATE)
    time_types = vcat(:maint, fill(:deg, n_deg), fill(:maint, n_maint))

    # Extract model parameters and data
    ρ = nlwienerard.efficiencies
    deg = Vector{Float64}(undef, n_deg)

    # Sort time and track types
    ordered_time_index = sortperm(time)
    ordered_time = time[ordered_time_index]
    ordered_time_types = time_types[ordered_time_index]

    # Find positions of maintenance and degradation dates in sorted array
    maint_index = findall(==(Symbol(:maint)), ordered_time_types)
    deg_index = findall(==(Symbol(:deg)), ordered_time_types)

    # Generate Wiener process values at all time points
    Y = nl_wiener_rand(nlwienerard, ordered_time)

    # Adjust for maintenance effects
    for i in eachindex(maint_index[1:end-1])
        if typeof(nlwienerard) == NLWienerARD1
                a = Y[maint_index[i + 1]] - Y[maint_index[i]]
        else
                a = Y[maint_index[i + 1]]
        end
        ρi = ρ[maint[i, "TYPE"]]
        Y[maint_index[i+1]:end] .-= ρi*a
    end

    deg = Y[deg_index]

    return deg
end