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
function mw_rand(mward::MWARD, inspection_dates::Vector{Float64})

    # parameters
    μ, σ = mward.drift, mward.volatility

    #time increments
    Δt = diff(inspection_dates)

    #simulation of the increments of the process according to the distribution of the increments of a Wiener process with drift
    X = [Δtij > 0 ? rand(MvNormal(Δtij * μ, sqrt(Δtij) * σ)) : zeros(length(μ)) for Δtij in Δt]

    #to have a matrix with each column corresponding to a time step
    reducedX = reduce(hcat, vcat([zeros(length(μ))], X))

    return  cumsum(reducedX, dims=2)
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
        if typeof(mward) == MWARD1
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

mutable struct MvDegradationData
    maintenances::DataFrame
    degradations::DataFrame    
end

function rand!(mward::MWARD, degradationdata::MvDegradationData)
    
    values = rand(mward, degradationdata.degradations.DATE, degradationdata.maintenances)

    for i in 1:length(mward.drift)
        println("VALUE$i")
        degradationdata.degradations[!, "VALUE$i"] = values[i, :]
    end

    return degradationdata
end