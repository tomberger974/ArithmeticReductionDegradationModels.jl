"""
    Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `time`.
    The simulation is done by simulating the increments of the process according to the distribution of the increments of a Wiener process with drift and then by taking the cumulative sum of these increments.
    Works only if the different dimensions are observed at the same time points.
"""
function mw_rand(mvw::MvWienerAR, inspection_dates::Vector{Float64})

    # parameters
    μ, Σ = mvw.drift, mvw.volatility

    #time increments
    Δt = diff(inspection_dates)

    #simulation of the increments of the process according to the distribution of the increments of a Wiener process with drift
    X = [Δtij > 0 ? rand(MvNormal(Δtij * μ, sqrt(Δtij) * Σ)) : zeros(length(μ)) for Δtij in Δt]

    #to have a matrix with each column corresponding to a time step
    reducedX = reduce(hcat, vcat([zeros(length(μ))], X))

    return cumsum(reducedX, dims=length(mvw.drift))
end

"""
    Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `inspection_dates` and with adjustments for maintenance effects according to the parameters of `mward` and the dates and types of maintenances given in `maintenances`.
    The simulation is done by simulating firstly an unmaintained process with `mw_rand` and then by adjusting the values of the process after each maintenance according to the type of maintenance and the corresponding efficiency given in `mward`.
"""
function rand(mvw::MvWienerAR, inspection_dates::Vector{Float64}, maintenances::DataFrame)

    # Create time vector with tags to track element types
    n_deg = length(inspection_dates)
    n_maint = nrow(maintenances)
    time = vcat(0., inspection_dates, maintenances.DATE)
    time_types = vcat(:maint, fill(:deg, n_deg), fill(:maint, n_maint))

    # Extract model parameters and data
    ρ = mvw.efficiencies
    r = length(mvw.drift)
    deg = Matrix{Float64}(undef, r, n_deg)
    indicators = unique(k[1] for k in keys(ρ))

    # Sort time and track types
    ordered_time_index = sortperm(time)
    ordered_time = time[ordered_time_index]
    ordered_time_types = time_types[ordered_time_index]

    # Find positions of maintenance and degradation dates in sorted array
    maint_index = findall(==(Symbol(:maint)), ordered_time_types)
    deg_index = findall(==(Symbol(:deg)), ordered_time_types)

    # Generate Wiener process values at all time points
    Y = mw_rand(mvw, ordered_time)

    # Adjust for maintenance effects
    for i in eachindex(maint_index[1:end-1])
        for j in eachindex(indicators)
            if typeof(ρ[(indicators[j], maintenances[i, "TYPE"])]) == ARD1
                a = Y[j, maint_index[i + 1]] - Y[j, maint_index[i]]
            else
                a = Y[j, maint_index[i + 1]]
            end
            ρui = ρ[(indicators[j], maintenances[i, "TYPE"])].value
            Y[j, maint_index[i+1]:end] .-= ρui * a
        end
    end

    # Extract degradation values at inspection dates to delete information related to maintenances
    deg = Y[:, deg_index]

    return deg
end

""" Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `inspection_dates` and with adjustments for maintenance effects according to the parameters of `mward` and the dates and types of maintenances given in `maintenances`.
    The simulation is done by simulating firstly an unmaintained process with `mw_rand` and then by adjusting the values of the process after each maintenance according to the type of maintenance and the corresponding efficiency given in `mward`.
"""
function rand!(mvw::MvWienerAR, degradationdata::DegradationData)
    
    values = rand(mvw, sort(unique(degradationdata.degradations.DATE)), degradationdata.maintenances)

    degradationdata.degradations.VALUE = reshape(values, :)

    return degradationdata
end

"""    delete_observations!(degradationdata::DegradationData)
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
    count_NB_MAINTENANCES(maintenance_dates::Vector{Float64}, inspection_dates::Vector{Float64})
    
    Return a vector of the same length as `inspection_dates` with the number of maintenances before each inspection date.
    This function is used to fill the column NB_MAINTENANCES of the DataFrame degradations of a DegradationData instance.
"""
function count_NB_MAINTENANCES(maintenance_dates::Vector{Float64}, inspection_dates::Vector{Float64})
    n = length(inspection_dates)
    NB_MAINTENANCES = zeros(Int, n)
    
    for i in 1:n
        NB_MAINTENANCES[i] = sum(maintenance_dates .< inspection_dates[i])
    end

    return NB_MAINTENANCES
end

"""
    DegradationData(mvw::MvWienerAR; K = 3, N_i = 5, Δt = 1.)
    
    Return a DegradationData instance with simulated degradation data according to the multidimensional Wiener process with drift and ARD maintenance model `mvw`.
    The simulation is done on a time grid of length `K` with `N_i` time steps between each maintenance, and a time increment of `Δt`.
"""
function DegradationData(mvw::MvWienerAR; K = 3, N_i = 5, Δt = 1., τ_types = rand(Set(k[2] for k in keys(mvw.efficiencies)), K), indicators = [Symbol("ind", i) for i in 1:length(mvw.drift)], deletion::Bool=false)
    r = length(mvw.drift)
    τ = convert(Vector{Float64}, [j*Δt*N_i for j in 1:K])
    T = convert(Float64, τ[end] + Δt*N_i)
    # Maintenances instance
    maintenances = DataFrame(DATE = τ, TYPE = τ_types)
    # Degradations with a single column VALUE and a column INDICATOR to precise what indicator does the value refers to
    degradations = DataFrame(DATE = Vector{Float64}([]), VALUE = Vector{Float64}([]), TYPE = Vector{Symbol}([]), NB_MAINTENANCES = Vector{Int64}([]))
    truc = vcat([0.], τ, T)
    for i in eachindex(truc)[1:end-1]
        row_after = DataFrame(DATE = [truc[i] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
        rows_between = DataFrame(DATE = vcat(fill(truc[i]+1:Δt:truc[i+1]-1, r)...), VALUE = Vector{Float64}(undef, (N_i - 1)*r), TYPE = vcat([[indicators[i] for _ in truc[i]+1:Δt:truc[i+1]-1] for i in 1:r]...), NB_MAINTENANCES = [i-1 for _ in (N_i - 1)*r])
        row_before = DataFrame(DATE = [truc[i+1] for _ in indicators], VALUE = Vector{Float64}(undef, r), TYPE = indicators, NB_MAINTENANCES = [i-1 for _ in r])
        degradations = vcat(degradations, row_after, rows_between, row_before)
    end
    # Degradation instance
    degradationdata = DegradationData(maintenances, degradations[r+1:end-r, :])
    # Deletion of a certain number of observations for a realistic case study
    # filter!(row -> row.DATE ∉ degradationdata.maintenances.DATE, degradationdata.degradations)
    # Deletion of some observations
    if deletion
        delete_observations!(degradationdata)
    end

    return degradationdata
end