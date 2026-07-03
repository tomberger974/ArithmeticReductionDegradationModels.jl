using Pkg
Pkg.activate(".")

using Revise

using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

using DataFrames
using Distributions
using LinearAlgebra

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

μ = [1., 1.]
Σ = diagm(ones(2))
ρ = Dict((:ind1, :M) => ARD.Efficiency(.5, ARD.ARD1()), (:ind2, :M) => ARD.Efficiency(.5, ARD.ARDinf()))
mvw = ARD.MvWienerAR(μ, Σ, ρ)
ARD.MvWienerAR()
keys(ρ)
Set(k[2] for k in keys(ρ))
rand(union(Set(k[2] for k in keys(ρ)), Set([:P])))

degradationdata = DegradationData(mvw)
degradations = degradationdata.degradations
maintenances = degradationdata.maintenances

ARD.mw_rand(mvw, sort(unique(degradationdata.degradations.DATE)))
ARD.rand!(mvw, degradationdata)

function truc(mvdegradationdata::MvDegradationData)
    deg = mvdegradationdata.degradations
    indicators = unique(deg.TYPE)
    n = nrow(filter(row -> row.TYPE ==:ind1, deg))

    indicators = unique(deg.TYPE)

    nb_inspections = rand(0:n)
    if nb_inspections == [0.]
        return mvdegradationdata
    end

    inspections = sample(1:n, nb_inspections; replace=(n == 1 ? true : false))

    for inspection in inspections
        for indicator in indicators
            rows = findall(row -> row.TYPE == indicator, deg)
            if isempty(rows)
                continue
            end

            row = rows[inspection]
            deg[row, :VALUE] = NaN
        end
    end

    return mvdegradationdata
end

A = [1, 3, 2]
sort(A)

rand!(mvw, degradationdata)
ARD.time_subdivision(mvdegradationdata)