module ArithmeticReductionDegradationModels

using DataFrames, Distributions, Random, CairoMakie, Optim, StatsBase, LinearAlgebra

export WienerARD, WienerARD∞, WienerARD1, params, params!, mean, var, std, quantile
export DegradationData, maintenances, degradations, maintenances!, degradations!, degradationsANDmaintenances!, infos, infos!
export rand, rand!
export MaintenancesPlot, DegradationsPlot, QuantilesPlot, plot!
export loglikelihood, fit_mle, fit_mle!, confint

include("dataset.jl")

include("unidim_models.jl")
include("unidim_tools_laurent.jl")
include("unidim_simulation.jl")
include("unidim_mle_ardinfinity.jl")
include("unidim_mle_ard1.jl")

include("nonlinear_unidim_models.jl")
include("nonlinear_unidim_simulation.jl")
include("nonlinear_unidim_mle_ardinfinity.jl")
include("nonlinear_unidim_mle_ard1.jl")

include("multidim_models.jl")
include("multidim_tools.jl")
include("multidim_simulation_dataset.jl")
include("multidim_mle.jl")

include("plot.jl")

end