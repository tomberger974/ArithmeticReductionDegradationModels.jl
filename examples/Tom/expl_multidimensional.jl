using Pkg
Pkg.activate(".")

using Revise

using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

using DataFrames
using Distributions
using LinearAlgebra

mward1 = ARD.MWARD1([1., 1.], diagm([1., 1.]), Dict("P" => [0.5, 0.2]))
mward∞ = ARD.MWARD∞([1., 1.], diagm([1., 1.]), Dict("P" => [0.5, 0.2]))
μ = mward1.drift
Σ = mward1.volatility
ρ = mward1.efficiencies

#degradationdata parameters
k = 3 #number of maintenances
nj = 5 #number of observations between each maintenance
τ = [i*nj for i in 1:k] #maintenance dates
τ_types = ["P" for i in 1:k] #maintenance types
T = τ[end] + nj #maximal time of observation
Δt = 1. #time increment between 2 observations

#create the degradationdata instance
degradationdata = DegradationData(τ, Δt, T, τ_types)
filter!(row -> row.TYPE in ["Between"], degradationdata.degradations)
maint = degradationdata.maintenances
deg = degradationdata.degradations


using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, deg.DATE, traj[1, :])
display(fig)


traj = ARD.rand(mward1, deg.DATE, maint)