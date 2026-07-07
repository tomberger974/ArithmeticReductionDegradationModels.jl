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

μ = [1., 1., 1.]
Σ = diagm(ones(3))
ρ = Dict((:ind1, :M) => ARD.Efficiency(.5, ARD.ARD1()), (:ind2, :M) => ARD.Efficiency(.5, ARD.ARDinf()), (:ind1, :C) => ARD.Efficiency(.4, ARD.ARD1()), (:ind2, :C) => ARD.Efficiency(.6, ARD.ARDinf()), (:ind2, :M) => ARD.Efficiency(.5, ARD.ARDinf()), (:ind3, :C) => ARD.Efficiency(.4, ARD.ARD1()), (:ind3, :C) => ARD.Efficiency(.6, ARD.ARDinf()))
mvw = ARD.MvWienerAR(μ, Σ, ρ)


degradationdata = ARD.DegradationData(mvw; deletion=false)
degradations = degradationdata.degradations
println(degradations)
maintenances = degradationdata.maintenances
ARD.rand!(mvw, degradationdata)


DataFrame(TRUC = 1.)

unique(filter(row -> row.NB_MAINTENANCES == 0, degradations).DATE)

ARD.count_inspections(degradationdata)