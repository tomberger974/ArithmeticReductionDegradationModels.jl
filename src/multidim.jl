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

"""
    count_inspections(degradationdata::DegradationData)

TBW
"""
function count_inspections(degradationdata::DegradationData)
    K = degradationdata.degradations[end, "NB_MAINTENANCES"]

    N_inspec = Vector{Int64}(undef, K + 1)

    for i in 0:K
        N_inspec[i+1] = length(unique(filter(row -> row.NB_MAINTENANCES == i, degradationdata.degradations).DATE))
    end

    return N_inspec
end


"""
    time_subdivision(deg::DataFrame, maint::DataFrame, i::Int64)

    Returns a serie of time increments for a given Inter-Maintenance-Interval (IMI), 
    going from the date of the maintenance action occuring at the start of the IMI 
    to the date of the maintenance action occuring at the end of the IMI.
    To be combined with time_subdivisions that gather each subdivision in a single dictionary.
"""
function time_subdivision(deg::DataFrame, maint::DataFrame, i::Int64)
    
    # find the inspection dates for which exactly i maintenance actions already occured
    maintenanceless_dates = sort(unique(filter(row -> row.NB_MAINTENANCES == i, deg).DATE))

    # find the maintenance dates occuring just before and just after maintenanceless_dates
    if i == 0
        maint_imoins1 = 0.
        maint_i = maint[i+1, "DATE"]
    elseif i == nrow(maint)
        maint_imoins1 = maint[i, "DATE"]
        maint_i = max(deg.DATE...)
    else
        maint_imoins1 = maint[i, "DATE"]
        maint_i = maint[i+1, "DATE"]
    end

    return diff(vcat(maint_imoins1, maintenanceless_dates, maint_i))
end


"""
    time_subdivisions(degradationdata::DegradationData)

    Returns a dictionary of Int64 going from 0 to the number of maintenance actions in degradationdata.
    Each entry i correspond to an Inter-Maintenance-Interval for which the dictionary returns 
    a serie of time increments going from the date of the maintenance action occuring at the start of the IMI 
    to the date of the maintenance action occuring at the end of the IMI.
"""
function time_subdivisions(deg::DataFrame, maint::DataFrame)
    return Dict(i => time_subdivision(deg, maint, i) for i in 0:nrow(maint))
end


function jump_matrix(degradationdata::DegradationData, mvw::MvWienerAR)

    # Data parameters
    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    # Extract model efficiencies and indicators
    ρ = mvw.efficiencies
    indicators = unique(key[1] for key in keys(ρ))

    # Respectively number of maintenance actions and indicators
    K = nrow(degradationdata.maintenances)

    # Provide a time subdivision for the matrices computation
    subdivision = time_subdivisions(deg, maint)
    subdivision_cumulative_length = cumsum([length(subdivision[i]) for i in 0:K])
    subdivision_total_length = subdivision_cumulative_length[K]

    # Vector containing the different blocks constituting the final matrix
    blocks = Vector{Matrix{Float64}}(undef, length(indicators))


    for p in eachindex(indicators)
        averaginginsertions = Vector{AveragingInsertion{Float64}}(undef, K)
        for i in 1:K
            ρuip = ρ[(indicators[p], maint[i, "TYPE"])]

            if ρuip.model isa ARD1
                averaginginsertions[i] = AveragingInsertion{Float64}(
                    subdivision_cumulative_length[i] + i,
                    i != 1 ? subdivision_cumulative_length[i-1] + i - 1 : 1, subdivision_cumulative_length[i] + i - 1,
                    ρuip.value)
            elseif ρuip.model isa ARDinf
                averaginginsertions[i] = AveragingInsertion{Float64}(
                    subdivision_cumulative_length[i] + i,
                    1, subdivision_cumulative_length[i] + i - 1,
                    ρuip.value)
            end
        end

        blocks[p] = build_block(subdivision_total_length + K, averaginginsertions)
    end

    return blocks
end





# An insertion defined by:
# - k          : position where the new row is inserted
# - first:last : columns of the non-zero coefficients in α
# - coefficient: common value of the non-zero coefficients
struct AveragingInsertion{T}
    k::Int
    first::Int
    last::Int
    coefficient::T
end


"""
    insert_row!(B, nrows, ins)

Apply one row-insertion transformation to the current block B.

The new row is:
    coefficient * sum(B[first:last, :], dims=1)

which corresponds to multiplying by the insertion row α.
"""
function insert_row!(B, nrows, ins::AveragingInsertion)

    newrow = ins.coefficient .* vec(sum(
        @view(B[ins.first:ins.last, :]),
        dims=1
    ))

    copyto!(
        @view(B[ins.k+1:nrows+1, :]),
        @view(B[ins.k:nrows, :])
    )

    B[ins.k, :] .= newrow

    return nrows + 1
end


"""
    build_block(n, insertions)

Construct the final block starting from an n×n identity matrix.

`insertions` is a vector of AveragingInsertion objects.
"""
function build_block(n::Int, insertions)

    # Final size is known beforehand
    B = zeros(Float64, n + length(insertions), n)

    # Initial identity
    B[1:n, :] .= Matrix(I, n, n)

    nrows = n

    for ins in insertions
        nrows = insert_row!(B, nrows, ins)
    end

    return B
end


function ARD_matrix(mvw::MvWienerAR, degradationdata::DegradationData)
    
end