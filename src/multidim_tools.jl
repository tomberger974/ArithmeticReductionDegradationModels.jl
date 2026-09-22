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


#MIGHT BE USELESS BUT NOT SURE ENOUGH TO DELETE
"""
    count_inspections(degradationdata::DegradationData)

Return a vector containing the number of inspections occurring in each Inter-Maintenance-Interval (IMI).
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
    time_subdivision(deg::DataFrame, maint::DataFrame, i::Int64; inspec_maint = false::Bool)

Returns a serie of time increments for a given Inter-Maintenance-Interval (IMI).
    
They go from the date of the maintenance action occuring at the start of the IMI to the date of the maintenance action occuring at the end of the IMI. 

If inspec_maint == true, inspections occuring at the same date as a maintenance action are not taken into account in the returned serie in order to avoid time increments of 0s which itself allow better computation performance.

To be combined with time_subdivisions that gather each subdivision in a single dictionary.
"""
function time_subdivision(deg::DataFrame, maint::DataFrame, i::Int64; inspec_maint = false::Bool)
    if inspec_maint
        # row.DATE ∉ maint.DATE in order to avoid time increments of 0s
        return sort(unique(filter(row -> row.NB_MAINTENANCES == i && row.DATE ∉ maint.DATE, deg).DATE))
    else
        return sort(unique(filter(row -> row.NB_MAINTENANCES == i, deg).DATE))
    end
end


"""
    time_subdivisions(degradationdata::DegradationData)

Returns a dictionary of Int64 going from 0 to the number of maintenance actions in degradationdata.

Each entry i correspond to an Inter-Maintenance-Interval for which the dictionary returns a serie of latent time increments going from the date of the maintenance action occuring at the start of the IMI to the date of the maintenance action occuring at the end of the IMI.
"""
function time_subdivisions(deg::DataFrame, maint::DataFrame; inspec_maint = false::Bool)
    return Dict(i => time_subdivision(deg, maint, i; inspec_maint = inspec_maint) for i in 0:nrow(maint))
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

The new row is coefficient * sum(B[first:last, :], dims=1).
"""
function insert_row!(B, nrows, ins::Insertion{Float64})

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

Construct a block starting from an n×n identity matrix.
"""
function build_block(n::Int, insertions::Vector{Insertion{Float64}})

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