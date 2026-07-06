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
        function verif(drift::Vector{Float64}, volatility::Matrix{Float64}, efficiencies::Efficiencies)
            s = length(drift)
            if size(volatility, 1) != s || size(volatility, 2) != s
                throw(ArgumentError("volatility must be square and match length(drift)=" * string(s)))
            end

            return true
        end

        verif(drift, volatility, efficiencies)
        new(drift, volatility, efficiencies)
    end
end

MvWienerAR(; drift::Vector{Float64}=zeros(Float64, 1), volatility::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 1, 1), efficiencies::Efficiencies=Dict((:ind1, :M) => Efficiency(.5, ARDinf()))) = MvWienerAR(drift, volatility, efficiencies)
# MvWienerAR(dim::Int64; efficiencies::Efficiencies) = MvWienerAR(zeros(dim), diagm(ones(dim)), efficiencies)

Base.show(io::IO, mvw::MvWienerAR) = print(io, "μ=", mvw.drift, ", Σ=", mvw.volatility, ", ρ=", mvw.efficiencies)

function time_subdivision(degradationdata::DegradationData)
    return sort(unique(vcat(degradationdata.degradations.DATE, degradationdata.maintenances.DATE)))
end

function count_inspections(degradationdata::DegradationData)
    K = degradationdata.degradations[end, "NB_MAINTENANCES"]

    N_inspec = Vector{Int64}(undef, K + 1)

    for i in 0:K
        N_inspec[i+1] = length(unique(filter(row -> row.NB_MAINTENANCES == i, degradationdata.degradations).DATE))
    end

    return N_inspec
end

function time_subdivision(degradationdata::DegradationData)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    inspection_dates = Vector{Float64}([])
    for i in 1:nrow(maint)
        push!(inspection_dates, unique(filter(row -> row.NB_MAINTENANCES == i, deg).DATE))
    end
    
    return sort(vcat([0.], inspections, maint.DATE))
end


function jump_matrix(degradationdata::DegradationData, mvw::MvWienerAR)

    # Data parameters
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    K = nrow(degradationdata.maintenances)

    # Create time vector with tags to track element types
    time = vcat(0., sort(unique(deg.DATE)), maint.DATE)
    time_types = vcat(:maint, fill(:deg, K), fill(:maint, nrow(degradationdata.maintenances)))

    # Extract model parameters and data
    ρ = mvw.efficiencies
    indicators = unique(k[1] for k in keys(ρ))

    # Sort time and track types
    ordered_time_index = sortperm(time)
    ordered_time = time[ordered_time_index]
    ordered_time_types = time_types[ordered_time_index]

    # Find positions of maintenance and degradation dates in sorted array
    maint_index = findall(==(Symbol(:maint)), ordered_time_types)
    deg_index = findall(==(Symbol(:deg)), ordered_time_types)

    for indicator in indicators
        B = I

        for i in 1:K
            ρuip = ρ[(indicator, maint[i, "TYPE"])]
            Δt = unique(filter(row -> row.NB_MAINTENANCES == i, deg).DATE)
            if ρuip.model isa ARD1
                ARD1_jump = vcat([0. for _ in 1:2])
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