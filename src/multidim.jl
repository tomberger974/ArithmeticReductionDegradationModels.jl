"""
    MaintenanceModels

Abstract type for Maintenance Models
"""
abstract type MaintenanceModels end


"""
    ARD1


"""
struct ARD1 <: MaintenanceModels end
struct ARDinf <: MaintenanceModels end


"""
    Efficiency

Efficiency object containing an imperfect maintenance model and a an efficiency value. Need to be linked 
with a tuple (indicator, maintenance type) to have meaning (cf. Efficiencies)
"""
mutable struct Efficiency
    value::Float64
    model::MaintenanceModels
end


"""
    Efficiencies = Union{
    Dict{Tuple{Symbol, Symbol}, Efficiency}, 
    Dict{Tuple{Symbol, String}, Efficiency}, 
    Dict{Tuple{String, Symbol}, Efficiency}, 
    Dict{Tuple{String, String}, Efficiency}}

Alias for a Union of different Dict in order to be able to use both string and symbols for the indicators 
or the maintenance type of the model. Contains a Dict linking each couple of indicator and maintenance 
type to an efficiency object containing an imperfect maintenance model and a an efficiency value

Order of arguments is mandatory: Efficiencies[indicator, maintenance_type]
"""
const Efficiencies = Union{
    Dict{Tuple{Symbol, Symbol}, Efficiency}, 
    Dict{Tuple{Symbol, String}, Efficiency}, 
    Dict{Tuple{String, Symbol}, Efficiency}, 
    Dict{Tuple{String, String}, Efficiency}}


"""
    MvWienerAR <: DPIM

Abstract type for Degradation Process with Imperfect Maintenance
"""
abstract type DPIM end


"""
    MvWienerAR <: DPIM

Type for Multivariate Wiener degradation process with imperfect maintenance following a model of type AR
"""
mutable struct MvWienerAR <: DPIM
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

Return a vector containing the number of inspections occurring in each Inter-Maintenance-Interval (IMI).
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
    return sort(unique(filter(row -> row.NB_MAINTENANCES == i, deg).DATE))
end


"""
    time_subdivisions(degradationdata::DegradationData)

Returns a dictionary of Int64 going from 0 to the number of maintenance actions in degradationdata.
Each entry i correspond to an Inter-Maintenance-Interval for which the dictionary returns a serie
of latent time increments going from the date of the maintenance action occuring at the start of the IMI
to the date of the maintenance action occuring at the end of the IMI.
"""
function time_subdivisions(deg::DataFrame, maint::DataFrame)
    return Dict(i => time_subdivision(deg, maint, i) for i in 0:nrow(maint))
end


"""
    ARDmatrix(degradationdata::DegradationData, mvw::MvWienerAR)

Returns a Dict of matrices associating each indicator to a matrix reconstructing the different 
virtual jumps.
"""
function ARDmatrix(degradationdata::DegradationData, mvw::MvWienerAR)

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
    subdivision_cumulative_length = cumsum([length(subdivision[i]) + 1 for i in 0:K])

    # Vector containing the different blocks constituting the final matrix
    blocks = Dict(ind => Matrix{Float64}(undef, subdivision_cumulative_length[end] + K, subdivision_cumulative_length[end]) for ind in indicators)

    for p in eachindex(indicators)
        virtualjumps = Vector{Insertion{Float64}}(undef, K)
        for i in 1:K
            ρuip = ρ[(indicators[p], maint[i, "TYPE"])]

            if ρuip.model isa ARD1
                virtualjumps[i] = Insertion{Float64}(
                    subdivision_cumulative_length[i] + i,
                    i != 1 ? subdivision_cumulative_length[i-1] + i : 1, subdivision_cumulative_length[i] + i - 1,
                    ρuip.value)
                    println(i != 1 ? subdivision_cumulative_length[i-1] + i - 1 : 1)
            elseif ρuip.model isa ARDinf
                virtualjumps[i] = Insertion{Float64}(
                    subdivision_cumulative_length[i] + i,
                    1, subdivision_cumulative_length[i] + i - 1,
                    ρuip.value)
            end
        end

        blocks[indicators[p]] = build_block(subdivision_cumulative_length[K] + K, virtualjumps)
    end

    return blocks
end





# An insertion defined by:
# - k          : position where the new row is inserted
# - first:last : columns of the non-zero coefficients in α
# - coefficient: common value of the non-zero coefficients
struct Insertion{T}
    k::Int
    first::Int
    last::Int
    coefficient::T
end


"""
    insert_row!(B, nrows, ins)

Apply one row-insertion transformation to the current block B.

The new row is coefficient * sum(B[first:last, :], dims=1) which corresponds to multiplying by the insertion row α.
"""
function insert_row!(B, nrows, ins::Insertion)

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

`insertions` is a vector of Insertion objects.
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


function observation_matrix(degradationdata::DegradationData, mvw::MvWienerAR)

    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    # Extract model efficiencies and indicators
    ρ = mvw.efficiencies
    indicators = unique(key[1] for key in keys(ρ))

    # Respectively number of maintenance actions and indicators
    K = nrow(maint)

    # Provide a time subdivision for the matrices computation
    subdivision = time_subdivisions(deg, maint)

    # Dict associating an indicator 
    all_indices = Dict(ind => [] for ind in indicators)

    for ind in indicators
        subdivision_p = time_subdivisions(filter(row -> row.TYPE == ind, deg), maint)

        count = 0

        summation_indices = Vector{Int64}(undef, sum([length(subdivision_p[i]) for i in 0:K]))
        index = 0

        for i in 0:K
            for truc in subdivision[i]
                count += 1
                if truc in subdivision_p[i]
                    index += 1
                    summation_indices[index] = count
                end
                count += 2
            end
        end
        
        all_indices[ind] = summation_indices
    end

    return all_indices

end

function combine_matrices(degradationdata::DegradationData, mvw::MvWienerAR)

    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    # Extract model efficiencies and indicators
    ρ = mvw.efficiencies
    indicators = unique(key[1] for key in keys(ρ))

    A = observation_matrix(degradationdata, mvw)
    B = ARDmatrix(degradationdata, mvw)
    C = Dict(ind => Matrix{Float64}(undef, length(A[ind]), size(B[ind], 2)) for ind in indicators)


    for ind in indicators
        Ap = vcat(1, A[ind])
        Bp = B[ind]

        for i in 1:length(Ap)-1
            C[ind][i, :] = sum([Bp[j, :] for j in Ap[i]:Ap[i+1]])
        end
    end

    return C
end