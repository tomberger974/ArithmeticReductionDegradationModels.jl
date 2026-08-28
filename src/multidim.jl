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



"""
    ARDmatrix(degradationdata::DegradationData, mvw::MvWienerAR)

Returns a Dict of matrices associating each indicator to a matrix B reconstructing the different latent jumps.
"""
function ARDmatrix(degradationdata::DegradationData, mvw::MvWienerAR)

    # Data parameters
    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    # Extract model efficiencies and indicators
    ρ = mvw.efficiencies
    indicators = collect(keys(mvw.drift))

    # Respectively number of maintenance actions and indicators
    K = nrow(degradationdata.maintenances)

    # Provide a time subdivision for the matrices computation
    subdivision = time_subdivisions(deg, maint; inspec_maint = true)
    subdivision_cumulative_length = cumsum([length(subdivision[i]) + (i != K ? 1 : 0) for i in 0:K])

    # Vector containing the different blocks constituting the final matrix
    blocks = Dict(ind => Matrix{Float64}(undef, subdivision_cumulative_length[end] + K, subdivision_cumulative_length[end]) for ind in indicators)

    for p in eachindex(indicators)
        latentjumps = Vector{Insertion{Float64}}(undef, K)
        for i in 1:K
            ρuip = ρ[(indicators[p], maint[i, "TYPE"])]

            if ρuip.model isa ARD1
                latentjumps[i] = Insertion{Float64}(
                    subdivision_cumulative_length[i] + i,
                    i != 1 ? subdivision_cumulative_length[i-1] + i : 1, subdivision_cumulative_length[i] + i - 1,
                    ρuip.value)
            elseif ρuip.model isa ARDinf
                latentjumps[i] = Insertion{Float64}(
                    subdivision_cumulative_length[i] + i,
                    1, subdivision_cumulative_length[i] + i - 1,
                    ρuip.value)
            end
        end

        blocks[indicators[p]] = build_block(subdivision_cumulative_length[end], latentjumps)
    end

    return blocks
end



"""
    observation_matrix(degradationdata::DegradationData, mvw::MvWienerAR)

Returns a dictionary linking each indicator to a serie of indices that allow to constitute the matrix A that constucts the observed VARIABLES out of the latent VARIABLES.

For instance if it returns [1, 5] for the indicator :ind, it means that we have a total of two observations for this indicator and that its first observed increment is equal to the first latent increment and that its second observed variable is equal to the sum of the 2nd, 3rd, 4th and 5th latent variables.
"""
function observation_matrix(degradationdata::DegradationData, mvw::MvWienerAR)

    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    # Extract model indicators
    indicators = collect(keys(mvw.drift))

    # Respectively number of maintenance actions and indicators
    K = nrow(maint)

    # Provide a time subdivisions for the matrices computation
    subdivisions = time_subdivisions(deg, maint; inspec_maint = false)

    # Dict associating an indicator 
    all_indices = Dict(ind => [] for ind in indicators)

    for ind in indicators
        subdivisions_p = time_subdivisions(filter(row -> row.TYPE == ind, deg), maint; inspec_maint = false)

        count = 0

        summation_indices = Vector{Int64}(undef, sum([length(subdivisions_p[i]) for i in 0:K]))
        index = 0

        for i in 0:K

            if isempty(subdivisions_p)
                continue
            elseif (i == 0 ? false : subdivisions[i][1] == maint[i,"DATE"])
                nothing
            else
                count += 1
            end

            for date in subdivisions[i]
                if date in subdivisions_p[i]
                    index += 1
                    summation_indices[index] = count
                end

                count += 1
            end

            if (i == K ? false : subdivisions[i][end] == maint[i+1,"DATE"])
                nothing
            else
                count += 1
            end
        end
        
        all_indices[ind] = summation_indices
    end

    return all_indices

end

"""
    combine_matrices(degradationdata::DegradationData, mvw::MvWienerAR)

Returns a dictionary linking each indicator to a matrix that construct the observed VARIABLES out of the latent INCREMENTS (pretty much the AB matrix).
"""
function combine_matrices(degradationdata::DegradationData, mvw::MvWienerAR)

    # Extract indicators
    indicators = collect(keys(mvw.drift))

    # Define the needed objects
    A = observation_matrix(degradationdata, mvw)
    B = ARDmatrix(degradationdata, mvw)
    C = Dict(ind => Matrix{Float64}(undef, length(A[ind]), size(B[ind], 2)) for ind in indicators)

    for ind in indicators
        Ap = A[ind]
        Bp = B[ind]

        for i in eachindex(Ap)
            C[ind][i, :] = sum([Bp[j, :] for j in (i != 1 ? (Ap[i-1]+1:Ap[i]) : 1:Ap[i])])
        end
    end

    return C
end

function observed_correlation_matrix(degradationdata::DegradationData, mvw::MvWienerAR)

    T = sort(unique(vcat(degradationdata.degradations.DATE, degradationdata.maintenances.DATE)))
    DT = Diagonal(T)

    indicators = collect(keys(mvw.drift))

    AB = combine_matrices(degradationdata, mvw)

    Σ_O = [mvw.covariances[(ind1, ind2)] * AB[ind1] * DT * transpose(AB[ind2]) for ind1 in indicators, ind2 in indicators]
    Σ_O = reduce(vcat, [reduce(hcat, Σ_O[i, :]) for i in axes(Σ_O, 1)])

    return Symmetric(Σ_O)
end

function observed_drift(degradationdata::DegradationData, mvw::MvWienerAR)

    indicators = collect(keys(mvw.drift))

    AB = combine_matrices(degradationdata, mvw)

    T = sort(unique(vcat(degradationdata.degradations.DATE, degradationdata.maintenances.DATE)))

    μ_O = reduce(vcat, [mvw.drift[ind] * AB[ind] * T for ind in indicators])

    return μ_O
end

function loglikelihood(degradationdata::DegradationData, mvw::MvWienerAR)
    μ_O = observed_drift(degradationdata, mvw)
    Σ_O = observed_correlation_matrix(degradationdata, mvw)

    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    indicators = collect(keys(mvw.drift))

    # Observed increments
    Y = reduce(vcat, [diff(vcat(0., sort(filter(row -> row.TYPE == ind, deg).VALUE))) for ind in indicators])

    return logpdf(MvNormal(μ_O, Σ_O), Y)
end

function fit_mle(degradationdata::DegradationData, mvw::MvWienerAR)
    fitted_mvw = deepcopy(mvw)
    return fit_mle!(degradationdata, fitted_mvw)
end

fit_mle(mvw::MvWienerAR, degradationdata::DegradationData) = fit_mle(degradationdata, mvw)

function fit_mle!(degradationdata::DegradationData, mvw::MvWienerAR)
    indicators = sort!(collect(keys(mvw.drift)), by=string)
    covariance_keys = [(ind1, ind2) for (i, ind1) in enumerate(indicators) for (j, ind2) in enumerate(indicators) if i <= j]
    efficiency_keys = sort!(collect(keys(mvw.efficiencies)), by=key -> (string(key[1]), string(key[2])))

    initial = vcat(
        [mvw.drift[ind] for ind in indicators],
        [mvw.covariances[key] for key in covariance_keys],
        [mvw.efficiencies[key].value for key in efficiency_keys]
    )

    drift_end = length(indicators)
    covariance_end = drift_end + length(covariance_keys)

    function update_model!(mvw, parameters)
        for (index, indicator) in enumerate(indicators)
            mvw.drift[indicator] = parameters[index]
        end
        for (index, (ind1, ind2)) in enumerate(covariance_keys)
            value = parameters[drift_end + index]
            mvw.covariances[(ind1, ind2)] = value
            mvw.covariances[(ind2, ind1)] = value
        end
        for (index, key) in enumerate(efficiency_keys)
            mvw.efficiencies[key].value = parameters[covariance_end + index]
        end
    end

    function objective(parameters)
        if any(!isfinite, parameters)
            return Inf
        end

        update_model!(mvw, parameters)

        try
            value = -loglikelihood(degradationdata, mvw)
            return isfinite(value) ? value : Inf
        catch error
            return error isa PosDefException ? Inf : rethrow()
        end
    end

    result = optimize(objective, initial)
    update_model!(mvw, result.minimizer)
    return mvw
end

fit_mle!(mvw::MvWienerAR, degradationdata::DegradationData) = fit_mle!(degradationdata, mvw)