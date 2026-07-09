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
ρ = Dict((:ind1, :M) => ARD.Efficiency(.5, ARD.ARD1()), (:ind2, :M) => ARD.Efficiency(.5, ARD.ARDinf()), (:ind1, :C) => ARD.Efficiency(.4, ARD.ARD1()), (:ind2, :C) => ARD.Efficiency(.6, ARD.ARDinf()), (:ind2, :M) => ARD.Efficiency(.5, ARD.ARDinf()))
mvw = ARD.MvWienerAR(μ, Σ, ρ)


degradationdata = ARD.DegradationData(mvw; deletion=false, before=false, after = false)
deg = degradationdata.degradations
maint = degradationdata.maintenances
ARD.rand!(mvw, degradationdata)

unique(deg, [:DATE, :NB_MAINTENANCES])
println(diff(ARD.unsorted_time_subdivision(degradationdata)[1]))
findall(x -> x == 1, [1, 2, 1, 3])


maint_dates = vcat([0.], sort(degradationdata.maintenances.DATE), max(deg.DATE...))
f(i::Int64) = diff(vcat(maint_dates[i+1], sort(unique(filter(row -> row.NB_MAINTENANCES == i, deg).DATE)), maint_dates[i+2]))
Dict(i => f(i) for i in 0:nrow(maint))