"""
    mw_rand(mvw::MvWienerAR, inspection_dates::Vector{Float64})

Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `time`.

The simulation is done by simulating the increments of the process according to the distribution of the increments of a Wiener process with drift and then by taking the cumulative sum of these increments.

Works only if the different dimensions are observed at the same time points.
"""
function mw_rand(mvw::MvWienerAR, inspection_dates::Vector{Float64})

    # parameters
    indicators = collect(keys(mvw.drift))
    μ = [mvw.drift[indicator] for indicator in indicators]
    Σ = [mvw.covariances[(ind1, ind2)] for ind1 in indicators, ind2 in indicators]

    # time increments 
    Δt = diff(vcat([0.], inspection_dates))

    # simulation of the increments of the process according to the distribution of the increments of a Wiener process with drift
    # defined with a condition to avoid the case where an element of Δt is equal to 0, which would lead to the covariance matrix being semi definite and thus non invertible
    normal_sample = [Δtij > 0 ? rand(MvNormal(Δtij * μ, sqrt(Δtij) * Σ)) : zeros(length(μ)) for Δtij in Δt]

    # to have a matrix with each column corresponding to a time step
    reduced_normal_sample = reduce(hcat, vcat([zeros(length(μ))], normal_sample))
    X_vect = cumsum(reduced_normal_sample, dims=2)[ :, 2:end]

    X_dict = Dict(indicators[p] => X_vect[p, :] for p in eachindex(indicators))

    return X_dict
end

"""
    rand(mvw::MvWienerAR, inspection_dates::Vector{Float64}, maintenances::DataFrame)

Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `inspection_dates` and with adjustments for maintenance effects according to the parameters of `mward` and the dates and types of maintenances given in `maintenances`.

The simulation is done by simulating firstly an unmaintained process with `mw_rand` and then by adjusting the values of the process after each maintenance according to the type of maintenance and the corresponding efficiency given in `mward`.

Simulate also the degradation level value just before the maintenance action and just after the maintenance action.

Assumes that the intersection of maintenance dates and maintenance dates is empty.
"""
function rand(mvw::MvWienerAR, inspection_dates::Vector{Float64}, maintenances::DataFrame; delete_before_after = false)

    # Create time vector with tags to track element types
    n_maint = nrow(maintenances)
    time = vcat(0., setdiff(inspection_dates, maintenances.DATE), maintenances.DATE, maintenances.DATE)
    time_types = vcat(:after, fill(:deg, length(setdiff(inspection_dates, maintenances.DATE))), fill(:before, n_maint), fill(:after, n_maint))

    # Sort time and track types
    ordered_time_index = sortperm(time)
    ordered_time = time[ordered_time_index]
    ordered_time_types = time_types[ordered_time_index]

    # Find positions of maintenance and degradation dates in the sorted array ordered_time
    before_index = findall(==(Symbol(:before)), ordered_time_types)
    after_index = findall(==(Symbol(:after)), ordered_time_types)
    deg_index = findall(==(Symbol(:deg)), ordered_time_types)

    # Generate Wiener process values at all time points
    Y = mw_rand(mvw, ordered_time)

    # Extract model parameters and data
    ρ = mvw.efficiencies

    # deg = Matrix{Float64}(undef, r, n_deg)
    indicators = collect(keys(mvw.drift))

    # Adjust for maintenance effects
    for indicator in indicators
        for i in eachindex(before_index)[1:end]
            if typeof(ρ[(indicator, maintenances[i, "TYPE"])].model) == ARD1
                a = Y[indicator][before_index[i]] - Y[indicator][after_index[i]]
            else
                a = Y[indicator][before_index[i]]
            end
            ρui = ρ[(indicator, maintenances[i, "TYPE"])].value
            println("aft", Y[indicator][after_index[i]])
            println("bef", Y[indicator][before_index[i] - 1])
            Y[indicator][after_index[i + 1]:end] .-= ρui * a
        end
    end

    if delete_before_after
        for indicator in indicators
            Y[indicator] = Y[indicator][deg_index]
        end
    end

    return Y
end

"""
    rand!(mvw::MvWienerAR, degradationdata::DegradationData)

Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `inspection_dates` and with adjustments for maintenance effects according to the parameters of `mward` and the dates and types of maintenances given in `maintenances`.

The simulation is done by simulating firstly an unmaintained process with `mw_rand` and then by adjusting the values of the process after each maintenance according to the type of maintenance and the corresponding efficiency given in `mward`.
"""
function rand!(mvw::MvWienerAR, degradationdata::DegradationData)
    
    maintenances = degradationdata.maintenances
    degradations = degradationdata.degradations
    indicators = collect(keys(mvw.drift))
    r = length(indicators)

    inspection_dates = time_subdivisions(degradations, maintenances; inspec_maint = true)

    # Degradations with a column VALUE and a column INDICATOR to precise what indicator does the value refers to
    degradations2 = DataFrame(DATE = Vector{Float64}([]), VALUE = Vector{Float64}([]), TYPE = Vector{Symbol}([]), NB_MAINTENANCES = Vector{Int64}([]))
    
    maint_dates = vcat([0.], maintenances.DATE, degradations.DATE[end])
    for i in eachindex(maint_dates)[1:end-1]
        row_after = DataFrame(DATE = [maint_dates[i] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
        Ni = length(inspection_dates[i - 1])
        rows_between = DataFrame(DATE = vcat([inspection_dates[i - 1] for _ in indicators]...), VALUE = Vector{Float64}(undef, Ni * r), TYPE = reshape([indicator for _ in inspection_dates[i - 1], indicator in indicators], Ni * r), NB_MAINTENANCES = [i-1 for _ in 1:Ni * r])
        if i != eachindex(maint_dates)[end-1]
            row_before = DataFrame(DATE = [maint_dates[i+1] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
        else
            row_before = similar(rows_between, 0)
        end
        degradations2 = vcat(degradations2, row_after, rows_between, row_before)
    end

    values = rand(mvw, unique(degradationdata.degradations).DATE, degradationdata.maintenances)

    for indicator in indicators
        degradations2[degradations2.TYPE .== indicator, "VALUE"] .= values[indicator]
    end

    # degradationdata.degradations.VALUE = reshape(values, :)

    return degradations2
end

"""
    delete_observations!(degradationdata::DegradationData)

Randomly delete some observations in the degradations DataFrame of `degradationdata`.

This function is used to test the ability of the estimation methods to handle unsynchonized data.
"""
function delete_observations!(degradationdata::DegradationData)
    deg = degradationdata.degradations
    indicators = unique(deg.TYPE)
    n = nrow(filter(row -> row.TYPE ==indicators[1], deg))

    indicators = unique(deg.TYPE)

    nb_inspections = rand(0:n)
    if nb_inspections == [0.]
        return degradationdata
    end

    inspections = sample(1:n, nb_inspections; replace=(n == 1 ? true : false))

    for inspection in deg[inspections, :DATE]
        nb_indicators = rand(1:length(indicators)-1)
        deleted_indicators = sample(indicators, nb_indicators; replace=(nb_indicators == 1 ? true : false))

        for indicator in deleted_indicators
            rows = findall(row -> row.DATE == inspection && row.TYPE in deleted_indicators, eachrow(deg))
            if isempty(rows)
                continue
            end

            delete!(deg, rows)
        end
    end

    return degradationdata
end

"""
    DegradationData(mvw::MvWienerAR; K = 3, N_i = 5, Δt = 1., τ_types = rand(Set(k[2] for k in keys(mvw.efficiencies)), K), indicators = [Symbol("ind", i) for i in 1:length(mvw.drift)], before::Bool=false, after::Bool=false, deletion::Bool=false)
    
Return a DegradationData instance with simulated degradation data according to the multidimensional Wiener process with drift and ARD maintenance model `mvw`.

The simulation is done on a time grid of length `K` with `N_i` time steps between each maintenance, and a time increment of `Δt`.

May contain observations just before or just after the maintenance.
"""
function DegradationData(mvw::MvWienerAR; K = 3, N_i = 5, Δt = 1., τ_types = rand(Set(k[2] for k in keys(mvw.efficiencies)), K), indicators = collect(keys(mvw.drift)), before::Bool=false, after::Bool=false, deletion::Bool=false)
    r = length(mvw.drift)
    τ = convert(Vector{Float64}, [j*Δt*(N_i + 1) for j in 1:K])
    T = convert(Float64, τ[end] + Δt*N_i)

    # Maintenances with a column DATE and a column VALUE
    maintenances = DataFrame(DATE = τ, TYPE = τ_types)
    # Degradations with a column VALUE and a column INDICATOR to precise what indicator does the value refers to
    degradations = DataFrame(DATE = Vector{Float64}([]), VALUE = Vector{Float64}([]), TYPE = Vector{Symbol}([]), NB_MAINTENANCES = Vector{Int64}([]))
    
    # Depending on the values of before and after create a degradations object with observations just before or just after maintenance actions
    maint_dates = vcat([0.], τ, T + Δt)
    if !after & !before
        for i in eachindex(maint_dates)[1:end-1]
            rows_between = DataFrame(DATE = vcat(fill(maint_dates[i]+Δt:Δt:maint_dates[i+1]-Δt, r)...), VALUE = Vector{Float64}(undef, N_i*r), TYPE = vcat([[indicators[i] for _ in maint_dates[i]+1:Δt:maint_dates[i+1]-1] for i in 1:r]...), NB_MAINTENANCES = [i-1 for _ in N_i*r])
            degradations = vcat(degradations, rows_between)
        end
    elseif !after & before
        for i in eachindex(maint_dates)[1:end-1]
            rows_between = DataFrame(DATE = vcat(fill(maint_dates[i]+Δt:Δt:maint_dates[i+1]-Δt, r)...), VALUE = Vector{Float64}(undef, N_i*r), TYPE = vcat([[indicators[i] for _ in maint_dates[i]+1:Δt:maint_dates[i+1]-1] for i in 1:r]...), NB_MAINTENANCES = [i-1 for _ in N_i*r])
            if i != eachindex(maint_dates)[end-1]
                row_before = DataFrame(DATE = [maint_dates[i+1] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
            else
                row_before = similar(rows_between, 0)
            end
            degradations = vcat(degradations, rows_between, row_before)
        end
    elseif after & !before
        for i in eachindex(maint_dates)[1:end-1]
            row_after = DataFrame(DATE = [maint_dates[i] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
            rows_between = DataFrame(DATE = vcat(fill(maint_dates[i]+Δt:Δt:maint_dates[i+1]-Δt, r)...), VALUE = Vector{Float64}(undef, N_i*r), TYPE = vcat([[indicators[i] for _ in maint_dates[i]+1:Δt:maint_dates[i+1]-1] for i in 1:r]...), NB_MAINTENANCES = [i-1 for _ in N_i*r])
            degradations = vcat(degradations, row_after, rows_between)
        end
    else
        for i in eachindex(maint_dates)[1:end-1]
            row_after = DataFrame(DATE = [maint_dates[i] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
            rows_between = DataFrame(DATE = vcat(fill(maint_dates[i]+Δt:Δt:maint_dates[i+1]-Δt, r)...), VALUE = Vector{Float64}(undef, N_i*r), TYPE = vcat([[indicators[i] for _ in maint_dates[i]+1:Δt:maint_dates[i+1]-1] for i in 1:r]...), NB_MAINTENANCES = [i-1 for _ in N_i*r])
            if i != eachindex(maint_dates)[end-1]
                row_before = DataFrame(DATE = [maint_dates[i+1] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
            else
                row_before = similar(rows_between, 0)
            end
            degradations = vcat(degradations, row_after, rows_between, row_before)
        end
    end

    # In the after case we just want to delete observations corresponding to the date 0.
    if after
        delete!(degradations, 1:r)
    end

    # Degradation instance
    degradationdata = DegradationData(maintenances, degradations)

    # Deletion of Ssome observations
    if deletion
        delete_observations!(degradationdata)
    end

    return degradationdata
end