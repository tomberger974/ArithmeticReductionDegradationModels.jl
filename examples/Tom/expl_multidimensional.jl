using Pkg
Pkg.activate(".")

using Revise

using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

using DataFrames
using LinearAlgebra

indicators = [:ind1, :ind2]
maintenance_types = [:M, :C]
μ = Dict(ind => 1. for ind in indicators)
Σ = Dict((ind1, ind2) => (ind1 == ind2 ? 1. : 0.) for ind1 in indicators, ind2 in indicators)
ρ = Dict((ind, mt) => ARD.Efficiency(.5, ARD.ARDinf()) for ind in indicators, mt in maintenance_types)
mvw = ARD.MvWienerAR(drift = μ, covariances = Σ, efficiencies = ρ)

degradationdata = ARD.DegradationData(mvw; K = 4, N_i = 5, τ_types = [:M, :C, :C, :M], deletion=false, before=false, after = false)
deg = degradationdata.degradations
maint = degradationdata.maintenances
ARD.rand!(mvw, degradationdata)
# filter!(row -> row.NB_MAINTENANCES != 0 || row.TYPE != :ind1, degradationdata.degradations)

ARD.ARDmatrix(degradationdata, mvw)[:ind1]
ARD.ARDmatrix(degradationdata, mvw)[:ind2]
ARD.ARDmatrix(degradationdata, mvw)

ARD.observation_matrix(degradationdata, mvw)
ARD.observation_matrix(degradationdata, mvw)[:ind1]
ARD.observation_matrix(degradationdata, mvw)[:ind2]

ARD.combine_matrices(degradationdata, mvw)
ARD.combine_matrices(degradationdata, mvw)[:ind1]
ARD.combine_matrices(degradationdata, mvw)[:ind2]

ARD.observed_correlation_matrix(degradationdata, mvw)
inv(ARD.observed_correlation_matrix(degradationdata, mvw))
ARD.observed_drift(degradationdata, mvw)
ARD.loglikelihood(degradationdata, mvw)
fit_mle(degradationdata, mvw)


### Test for the insertion functions inside of matrices ###
ins1 = ARD.Insertion(3, 1, 2, .5)
ins2 = ARD.Insertion(5, 1, 3, .5)
ARD.build_block(4, [ins1, ins2])


### Test time_subdivision_tool ###
# Usual case
degradationdata = ARD.DegradationData(mvw; deletion=false, before=false, after = false)
ARD.time_subdivision(degradationdata.degradations, degradationdata.maintenances, 0)
ARD.time_subdivisions(degradationdata.degradations, degradationdata.maintenances; inspec_maint = false)
ARD.time_subdivisions(degradationdata.degradations, degradationdata.maintenances; inspec_maint = true)
# No observations before the first maintenance action
degradationdata = ARD.DegradationData(mvw; deletion=false, before=false, after = false)
filter!(row -> row.NB_MAINTENANCES != 0, degradationdata.degradations)
ARD.time_subdivisions(degradationdata.degradations, degradationdata.maintenances)



### Test for matrix builder ###
# without both before and after
degradationdata = ARD.DegradationData(mvw; K = 2, N_i = 1, τ_types = [:M, :C ], deletion=false, before = false, after = false)
ARD.rand!(mvw, degradationdata)
filter!(row -> row.NB_MAINTENANCES != 0 || row.TYPE != :ind1, degradationdata.degradations)

ARD.ARDmatrix(degradationdata, mvw)
ARD.ARDmatrix(degradationdata, mvw)[:ind1]
ARD.ARDmatrix(degradationdata, mvw)[:ind2]

ARD.observation_matrix(degradationdata, mvw)
ARD.observation_matrix(degradationdata, mvw)[:ind1]
ARD.observation_matrix(degradationdata, mvw)[:ind2]

ARD.combine_matrices(degradationdata, mvw)
ARD.combine_matrices(degradationdata, mvw)[:ind1]
ARD.combine_matrices(degradationdata, mvw)[:ind2]


# with before and after
degradationdata = ARD.DegradationData(mvw; K = 2, N_i = 1, τ_types = [:M, :C ], deletion=false, before=true, after = true)
ARD.rand!(mvw, degradationdata)
filter!(row -> row.NB_MAINTENANCES != 0 || row.TYPE != :ind1, degradationdata.degradations)

ARD.ARDmatrix(degradationdata, mvw)
ARD.ARDmatrix(degradationdata, mvw)[:ind1]
ARD.ARDmatrix(degradationdata, mvw)[:ind2]

ARD.observation_matrix(degradationdata, mvw)
ARD.observation_matrix(degradationdata, mvw)[:ind1]
ARD.observation_matrix(degradationdata, mvw)[:ind2]

ARD.combine_matrices(degradationdata, mvw)
ARD.combine_matrices(degradationdata, mvw)[:ind1]
ARD.combine_matrices(degradationdata, mvw)[:ind2]


# with before only
degradationdata = ARD.DegradationData(mvw; K = 2, N_i = 1, τ_types = [:M, :C ], deletion=false, before=true, after = false)
ARD.rand!(mvw, degradationdata)
filter!(row -> row.NB_MAINTENANCES != 0 || row.TYPE != :ind1, degradationdata.degradations)

ARD.ARDmatrix(degradationdata, mvw)
ARD.ARDmatrix(degradationdata, mvw)[:ind1]
ARD.ARDmatrix(degradationdata, mvw)[:ind2]

ARD.observation_matrix(degradationdata, mvw)
ARD.observation_matrix(degradationdata, mvw)[:ind1]
ARD.observation_matrix(degradationdata, mvw)[:ind2]

ARD.combine_matrices(degradationdata, mvw)
ARD.combine_matrices(degradationdata, mvw)[:ind1]
ARD.combine_matrices(degradationdata, mvw)[:ind2]


# with after only
degradationdata = ARD.DegradationData(mvw; K = 2, N_i = 1, τ_types = [:M, :C ], deletion=false, before=false, after = true)
ARD.rand!(mvw, degradationdata)
filter!(row -> row.NB_MAINTENANCES != 0 || row.TYPE != :ind1, degradationdata.degradations)

ARD.ARDmatrix(degradationdata, mvw)
ARD.ARDmatrix(degradationdata, mvw)[:ind1]
ARD.ARDmatrix(degradationdata, mvw)[:ind2]

ARD.observation_matrix(degradationdata, mvw)
ARD.observation_matrix(degradationdata, mvw)[:ind1]
ARD.observation_matrix(degradationdata, mvw)[:ind2]

ARD.combine_matrices(degradationdata, mvw)
ARD.combine_matrices(degradationdata, mvw)[:ind1]
ARD.combine_matrices(degradationdata, mvw)[:ind2]
