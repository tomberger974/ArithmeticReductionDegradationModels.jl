abstract type MaintenanceModels end
mutable struct ARD1 <: MaintenanceModels end
mutable struct ARDinf <: MaintenanceModels end

abstract type PDAMI end

mutable struct Efficiency
    value::Float64
    model::MaintenanceModels
end

const Efficiencies = Union{
    Dict{Tuple{Symbol, Symbol}, Efficiency}, 
    Dict{Tuple{Symbol, String}, Efficiency}, 
    Dict{Tuple{String, Symbol}, Efficiency}, 
    Dict{Tuple{String, String}, Efficiency}
}

mutable struct MvWienerAR <: PDAMI
    drift::Vector{Float64}
    volatility::Matrix{Float64}
    efficiencies::Efficiencies

    function MvWienerAR(drift::Vector{Float64}, volatility::Matrix{Float64}, efficiencies::Efficiencies)
        s = length(drift)
        if size(volatility, 1) != s || size(volatility, 2) != s
            throw(ArgumentError("volatility must be square and match length(drift)=" * string(s)))
        end
        new(drift, volatility, efficiencies)
    end
end

MvWienerAR(; drift::Vector{Float64}=zeros(Float64, 1), volatility::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 1, 1), efficiencies::Efficiencies=Dict((:ind1, :M) => Efficiency(.5, ARDinf()))) = MvWienerAR(drift, volatility, efficiencies)
MvWienerAR(dim::Int64; efficiencies::Efficiencies) = MvWienerAR(zeros(dim), diagm(ones(dim)), efficiencies)

Base.show(io::IO, mvw::MvWienerAR) = print(io, "μ=", mvw.drift, ", Σ=", mvw.volatility, ", ρ=", mvw.efficiencies)

function count_inspections(degradationdata::DegradationData)
    K = degradationdata.degradations[end, "NB_MAINTENANCES"]

    N_inspec = Vector{Int64}(undef, K + 1)

    for i in 0:K
        N_inspec[i+1] = length(unique(filter(row -> row.NB_MAINTENANCES == i, degradationdata.degradations).DATE))
    end

    return N_inspec
end

"""
    unsorted_time_subdivision(degradationdaa::DegradationData)
    Return two vectors: the first one merges all inspections and maintenances dates along with the 0 date
    and the second one indicates for each date if it is a maintenance or an inspection date with the convention that date 0 correspond to a maintenance.
"""
function unsorted_time_subdivision(degradationdata::DegradationData)
    deg = degradationdata.degradations
    maint_dates = vcat([0.], sort(degradationdata.maintenances.DATE), max(deg.DATE...))

    # vector with the time subdivision of interest
    time = Vector{Float64}([])
    # nb_maintenances[i] indicates how many maintenances occured at time[i]
    nb_maintenances = Vector{Int64}([])
    # inspec_or_maint[i] indicates if time[i] correspond to a maintenance date or to an inspection date
    inspec_or_maint = Vector{Symbol}([])
    
    for i in eachindex(maint_dates)[1:end-1]
        ith_deg = sort(unique(filter(row -> row.NB_MAINTENANCES == i-1, deg).DATE))
        time = vcat(time, maint_dates[i], ith_deg, maint_dates[i+1])
        nb_maintenances = vcat(nb_maintenances, [i for _ in 1:length(ith_deg)+2])
        inspec_or_maint = vcat(inspec_or_maint, :maint, fill(:deg, length(ith_deg)), :maint)
        truc = Dict()
    end

    return time, nb_maintenances, inspec_or_maint
end


I(n) = 
function jump_matrix(degradationdata::DegradationData, mvw::MvWienerAR)

    # Data parameters
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    K = nrow(degradationdata.maintenances)

    # Extract model parameters and data
    ρ = mvw.efficiencies
    indicators = unique(k[1] for k in keys(ρ))

    time, nb_maintenances, inspec_or_maint = unsorted_time_subdivision(degradationdata::DegradationData)

    f(i::Int64) = diff(vcat(maint_dates[i+1], sort(unique(filter(row -> row.NB_MAINTENANCES == i, deg).DATE)), maint_dates[i+2]))
    machin = Dict(i => f(i) for i in 0:nrow(maint))
    machin_length = Dict(i => length(machin(i)))
    machin_total_length = sum(machin_length(i))


    for indicator in indicators
        B = I

        for i in 0:K
            ρuip = ρ[(indicator, maint[i, "TYPE"])]

            if ρuip.model isa ARD1
                n = 5
                Iₙ = Matrix(I, n, n)

                r = [10, 20, 30, 40, 50]   # new row
                k = 3       
                B *= vcat
            end
        end
        B[1,1] = i
        # ...

        push!(blocks, B)
    end

    BD = BlockDiagonal(blocks)

    return jump_matrix

end

function machin(degrdationdata::DegradationData)
    deg = unique(deg, [:DATE, :NB_MAINTENANCES])
    
end