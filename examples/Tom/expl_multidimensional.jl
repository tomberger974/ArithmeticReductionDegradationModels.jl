using Pkg
Pkg.activate(".")

using Revise

using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

using DataFrames
using LinearAlgebra

μ = [1., 1.]
Σ = diagm(ones(2))
ρ = Dict((:ind1, :M) => ARD.Efficiency(.1, ARD.ARD1()), (:ind2, :M) => ARD.Efficiency(.1, ARD.ARDinf()), (:ind1, :C) => ARD.Efficiency(.1, ARD.ARD1()), (:ind2, :C) => ARD.Efficiency(.1, ARD.ARDinf()))
mvw = ARD.MvWienerAR(μ, Σ, ρ)


degradationdata = ARD.DegradationData(mvw; K = 2, N_i = 2, τ_types = [:M, :C ], deletion=false, before=false, after = false)
deg = degradationdata.degradations
maint = degradationdata.maintenances
ARD.rand!(mvw, degradationdata)
filter!(row -> row.NB_MAINTENANCES != 0 || row.TYPE != :ind1, degradationdata.degradations)

ARD.ARDmatrix(degradationdata, mvw)[:ind1]
ARD.ARDmatrix(degradationdata, mvw)[:ind2]
ARD.ARDmatrix(degradationdata, mvw)

ARD.observation_matrix(degradationdata, mvw)
ARD.observation_matrix(degradationdata, mvw)[:ind1]
ARD.observation_matrix(degradationdata, mvw)[:ind2]

ARD.combine_matrices(degradationdata, mvw)
ARD.combine_matrices(degradationdata, mvw)[:ind1]
ARD.combine_matrices(degradationdata, mvw)[:ind2]


#Test for the insertion functions inside of matrices
x = 1/2

ins1 = ARD.Insertion(
    3,      # insert after row 2
    1, 2,   # non-zero coefficients on columns 1:2
    x
)

y = 1/3

ins2 = ARD.Insertion(
    5,      # insert after row 4 of the current block
    1, 3,   # non-zero coefficients on columns 1:3
    y
)

ARD.build_block(4, [ins1, ins2])



### Test time_subdivision_tool ###
# Usual case
degradationdata = ARD.DegradationData(mvw; deletion=false, before=false, after = false)
ARD.time_subdivision(degradationdata.degradations, degradationdata.maintenances, 0)
ARD.time_subdivisions(degradationdata)
# No observations before the first maintenance action
degradationdata = ARD.DegradationData(mvw; deletion=false, before=false, after = false)
filter!(row -> row.NB_MAINTENANCES != 0, degradationdata.degradations)
ARD.time_subdivision(degradationdata.degradations, degradationdata.maintenances, 0)