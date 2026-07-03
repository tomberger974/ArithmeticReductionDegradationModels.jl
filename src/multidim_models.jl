abstract type MWARD end

mutable struct MWARD∞ <: MWARD
    drift::Vector{Float64}
    volatility::Matrix{Float64}
    efficiencies::Dict{Union{String, Symbol}, Vector{Float64}}
end

MWARD∞(; drift::Vector{Float64}=[0.], volatility::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 1, 1), efficiencies::Dict{Symbol, Float64}=Dict(:M => 0.)) = MWARD∞(drift, volatility, [efficiencies])
MWARD∞(dim::Int64;) = MWARD∞(zeros(dim), diagm(ones(dim)), [Dict(:M => 0.) for _ in 1:dim])

mutable struct MWARD1 <: MWARD
    drift::Vector{Float64}
    volatility::Matrix{Float64}
    efficiencies::Dict{Union{String, Symbol}, Vector{Float64}}
end

MWARD1(; drift::Vector{Float64}=[0.], volatility::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 1, 1), efficiencies::Dict{Symbol, Float64}=Dict(:M => 0.)) = MWARD1(drift, volatility, [efficiencies])
MWARD1(dim::Int64;) = MWARD1(zeros(dim), diagm(ones(dim)), [Dict(:M => 0.) for _ in 1:dim])

Base.show(io::IO, mward::MWARD) = print(io, "drift:", mward.drift, ", volatility:", mward.volatility, ", efficiencies:", mward.efficiencies)


"""
    Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `time`.
    The simulation is done by simulating the increments of the process according to the distribution of the increments of a Wiener process with drift and then by taking the cumulative sum of these increments.
    Works only if the different dimensions are observed at the same time points.
"""
function mw_rand(mvw::MvW, inspection_dates::Vector{Float64})

    # parameters
    μ, σ = mvw.drift, mvw.volatility

    #time increments
    Δt = diff(inspection_dates)

    #simulation of the increments of the process according to the distribution of the increments of a Wiener process with drift
    X = [Δtij > 0 ? rand(MvNormal(Δtij * μ, sqrt(Δtij) * σ)) : zeros(length(μ)) for Δtij in Δt]

    #to have a matrix with each column corresponding to a time step
    reducedX = reduce(hcat, vcat([zeros(length(μ))], X))

    return cumsum(reducedX, dims=length(mvw.drift))
end


"""
    Simulate a wiener process with drift according to the parameters of `mward` at the time points given in `inspection_dates` and with adjustments for maintenance effects according to the parameters of `mward` and the dates and types of maintenances given in `maintenances`.
    The simulation is done by simulating firstly an unmaintained process with `mw_rand` and then by adjusting the values of the process after each maintenance according to the type of maintenance and the corresponding efficiency given in `mward`.
"""
function rand(mward::MWARD, inspection_dates::Vector{Float64}, maintenances::DataFrame)

    # Create time vector with tags to track element types
    n_deg = length(inspection_dates)
    n_maint = nrow(maintenances)
    time = vcat(0., inspection_dates, maintenances.DATE)
    time_types = vcat(:maint, fill(:deg, n_deg), fill(:maint, n_maint))

    # Extract model parameters and data
    ρ = mward.efficiencies
    r = length(mward.drift)
    deg = Matrix{Float64}(undef, r, n_deg)

    # Sort time and track types
    ordered_time_index = sortperm(time)
    ordered_time = time[ordered_time_index]
    ordered_time_types = time_types[ordered_time_index]

    # Find positions of maintenance and degradation dates in sorted array
    maint_index = findall(==(Symbol(:maint)), ordered_time_types)
    deg_index = findall(==(Symbol(:deg)), ordered_time_types)

    # Generate Wiener process values at all time points
    Y = mw_rand(mward, ordered_time)

    # Adjust for maintenance effects
    for i in eachindex(maint_index[1:end-1])
        if typeof(mvw) == MWARD1
            a = Y[:, maint_index[i + 1]] - Y[:, maint_index[i]]
        else
            a = Y[:, maint_index[i + 1]]
        end
        ρi = ρ[maintenances[i, "TYPE"]]
        Y[:, maint_index[i+1]:end] .-= ρi .* a
    end

    # Extract degradation values at inspection dates to delete information related to maintenances
    deg = Y[:, deg_index]

    return deg
end

function rand(mvw::MvW, inspection_dates::Vector{Float64}, maintenances::DataFrame)

    # Create time vector with tags to track element types
    n_deg = length(inspection_dates)
    n_maint = nrow(maintenances)
    time = vcat(0., inspection_dates, maintenances.DATE)
    time_types = vcat(:maint, fill(:deg, n_deg), fill(:maint, n_maint))

    # Extract model parameters and data
    ρ = mvw.efficiencies
    r = length(mvw.drift)
    deg = Matrix{Float64}(undef, r, n_deg)

    # Sort time and track types
    ordered_time_index = sortperm(time)
    ordered_time = time[ordered_time_index]
    ordered_time_types = time_types[ordered_time_index]

    # Find positions of maintenance and degradation dates in sorted array
    maint_index = findall(==(Symbol(:maint)), ordered_time_types)
    deg_index = findall(==(Symbol(:deg)), ordered_time_types)

    # Generate Wiener process values at all time points
    Y = mw_rand(mvw, ordered_time)
    println(size(Y))

    # Adjust for maintenance effects
    for i in eachindex(maint_index[1:end-1])
        if typeof(ρ[maintenances[i, "TYPE"]].models) == ARD1
            a = Y[:, maint_index[i + 1]] - Y[:, maint_index[i]]
        else
            a = Y[:, maint_index[i + 1]]
        end
        ρi = ρ[maintenances[i, "TYPE"]].values
        Y[:, maint_index[i+1]:end] .-= ρi .* a
    end

    # Extract degradation values at inspection dates to delete information related to maintenances
    println(length(deg))
    deg = Y[:, deg_index]

    return deg
end

function rand!(mvw::MvW, degradationdata::DegradationData)
    
    values = rand(mvw, sort(unique(degradationdata.degradations.DATE)), degradationdata.maintenances)

    degradationdata.degradations.VALUE = reshape(values, :)

    return degradationdata
end


"""    delete_observations!(degradationdata::MvDegradationData)
    Randomly delete some observations in the degradations DataFrame of `degradationdata`.
    This function is used to test the ability of the estimation methods to handle unsynchonized data.
"""
function delete_observations!(mvdegradationdata::MvDegradationData)
    deg = mvdegradationdata.degradations
    indicators = unique(deg.TYPE)
    n = nrow(filter(row -> row.TYPE ==indicators[1], deg))

    indicators = unique(deg.TYPE)

    nb_inspections = rand(0:n)
    if nb_inspections == [0.]
        return mvdegradationdata
    end

    inspections = sample(1:n, nb_inspections; replace=(n == 1 ? true : false))

    for inspection in deg[inspections, :DATE]
        nb_indicators = rand(1:length(indicators))
        deleted_indicators = sample(indicators, nb_indicators; replace=(nb_indicators == 1 ? true : false))

        for indicator in deleted_indicators
            println(inspection)
            rows = findall(row -> row.DATE == inspection && row.TYPE in deleted_indicators, eachrow(deg))
            if isempty(rows)
                continue
            end

            delete!(deg, rows)
        end
    end

    return mvdegradationdata
end

function create_degradationdata(mvw::MvW; K = 3, N_i = 5, Δt = 1., τ = convert(Vector{Float64}, [j*Δt*N_i for j in 1:K]), τ_types = rand(keys(mvw.efficiencies), K),     T = convert(Float64, τ[end] + Δt*N_i))
    # Degradations with a single column VALUE and a column INDICATOR to precise what indicator does the value refers to
    degradations = DataFrame(DATE = vcat(1:Δt:T, 1:Δt:T), VALUE = Vector{Float64}(undef, 2 * Int64(T/Δt)), TYPE = vcat([:ind1 for _ in 1:Int64(T/Δt)], [:ind2 for _ in 1:Int64(T/Δt)]))
    # Maintenances instance
    maintenances = DataFrame(DATE = τ, TYPE = τ_types)
    # Degradation instance
    mvdegradationdata = MvDegradationData(maintenances, degradations)
    # Deletion of a certain number of observations for a realistic case study
    filter!(row -> row.DATE ∉ mvdegradationdata.maintenances.DATE, mvdegradationdata.degradations)
    # Simulation of an ARD process
    rand!(mvw, mvdegradationdata)
    # Deletion of some observations
    delete_observations!(mvdegradationdata)

    return DegradationData
end

function compute_N_inspec(mvdegradationdata::MvDegradationData)
    n = mvdegradationdata.degradations[end, "NB_MAINTENANCES"]

    N_inspec = Vector{Int64}(undef, n + 1)

    for i in 0:n
        N_inspec[i+1] = count(==(i), mvdegradationdata.degradations.NB_MAINTENANCES)
    end

    return N_inspec
end

abstract type ARD end
mutable struct ARDinf <: ARD end
mutable struct ARD1 <: ARD end

mutable struct efficiency
    models::Vector{ARD}
    values::Vector{Float64}
end

function efficiency(models::Vector{<:ARD}, values::Vector{Float64})
    efficiency(Vector{ARD}(models), values)
end

efficiency(; models::Vector{<:ARD}=[ARDinf()], values::Vector{Float64}=[0.]) = efficiency(models, values)
efficiency(dim::Int64; models::Vector{<:ARD}=fill(ARD1(), dim)) = efficiency(models, zeros(Float64, dim))

mutable struct MvW
    drift::Vector{Float64}
    volatility::Matrix{Float64}
    efficiencies::Dict{Union{String, Symbol}, efficiency}
end

MvW(; drift::Vector{Float64}=zeros(Float64, 1), volatility::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 1, 1), efficiencies::Dict{Symbol, efficiency}=Dict(:M => efficiency())) = MvW(drift, volatility, efficiencies)
MvW(dim::Int64; models::Vector{<:ARD}=fill(ARD1(), dim)) = MvW(zeros(dim), diagm(ones(dim)), Dict(:M => efficiency(dim; models)))

function time_subdivision(mvdegradationdata::MvDegradationData)
    return sort(unique(vcat(mvdegradationdata.degradations.DATE, mvdegradationdata.maintenances.DATE)))
end