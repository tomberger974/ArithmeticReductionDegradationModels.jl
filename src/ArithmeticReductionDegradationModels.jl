module ArithmeticReductionDegradationModels

using DataFrames, Distributions, Random, CairoMakie, Optim, StatsBase, LinearAlgebra

export WienerARD, WienerARD∞, WienerARD1, params, params!, mean, var, std, quantile
export DegradationData, maintenances, degradations, maintenances!, degradations!, degradationsANDmaintenances!, infos, infos!
export rand, rand!
export MaintenancesPlot, DegradationsPlot, QuantilesPlot, plot!
export loglikelihood, fit_mle, fit_mle!, confint

include("Tools.jl")
include("DataSet.jl")
include("Models.jl")
include("Simulation.jl")
include("Plot.jl")
include("mle_ardinfinity.jl")
include("mle_ard1.jl")

#include("multidim_models.jl")
include("multidim.jl")
include("multidim_data_simu.jl")

include("nl_simulation.jl")
include("nl_mle.jl")


end