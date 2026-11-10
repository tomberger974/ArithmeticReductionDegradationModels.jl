using Pkg

Pkg.activate(".")

using Revise

using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

using DataFrames
using LinearAlgebra
using Random
using CairoMakie

### PARAMETERS AND SETTINGS OF THE SIMULATIONS ###

# Names of the indicators and the maintenance types
indicators = [:ind1, :ind2]
maintenance_types = ["P"]

# Parameters of the multidimensional wiener process
μ = Dict(:ind1 => 1., :ind2 => 2.)
c11, c22 = .4, .2
c12 = .2 * (sqrt(c11) * sqrt(c22))
Σ = Dict((:ind1, :ind1) => c11, (:ind2, :ind2) => c22, (:ind1, :ind2) => c12, (:ind2, :ind1) => c12)

# Maintenance parameters for ARD models
ρ = Dict((:ind1, "P") => ARD.Efficiency(.5, ARD.ARDinf()), (:ind2, "P") => ARD.Efficiency(.5, ARD.ARDinf()))

# Instantiate the Multivariate Wiener ARD process
true_mvw = ARD.MvWienerAR(drift = μ, covariances = Σ, efficiencies = ρ)

# Instantiate the data
τ_types = ["P" for _ in 1:5]
degradationdata = ARD.DegradationData(true_mvw; nb_maint = 10, nb_inspec_per_IMI = 3, deletion=false, before=false, after = false)
deg = degradationdata.degradations
maint = degradationdata.maintenances

# Choose a total number of simulations N and a random seed for reproductability
N = 5000
nb_seed = 3
Random.seed!(nb_seed)



############## SIMULATIONS ##################
truc = ARD.tests(true_mvw, degradationdata, N)
#############################################



### PUT EVERYTHING IN VECTORS FOR THEM TO BE MANIPULATED MORE EASILY ###

# Drift estimations vectors
vec_drift1 = ARD.param.(truc, :drift; sel = :ind1)
vec_drift2 = ARD.param.(truc, :drift; sel = :ind2)

# Variance and covariance estimations vectors
vec_cov11 = ARD.param.(truc, :covariances, sel = (:ind1, :ind1))
vec_cov12 = ARD.param.(truc, :covariances, sel = (:ind1, :ind2))
vec_cov21 = ARD.param.(truc, :covariances, sel = (:ind2, :ind1))
vec_cov22 = ARD.param.(truc, :covariances, sel = (:ind2, :ind2))

# Efficiencies estimations vectors
vec_eff1P = getproperty.(ARD.param.(truc, :efficiencies, sel = (:ind1, "P")), :value)
vec_eff2P = getproperty.(ARD.param.(truc, :efficiencies, sel = (:ind2, "P")), :value)

sum()

### BOXPLOTS TO REPRESENT THE ESTIMATIONS

# Drift boxplots
fig = Figure()
    ax1 = Axis(
        fig[1,1],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = "drift 1"
        )
    ax2 = Axis(
        fig[1,2],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = "drift 2"
        )
    boxplot!(ax1, ones(length(vec_drift1)),vec_drift1)
    boxplot!(ax2, ones(length(vec_drift2)),vec_drift2)
    linkyaxes!(ax1, ax2)
    fig
    save("/Users/berger/Documents/latex/pics//multidim_simu1_drift_box.png", fig)



# Covariances boxplots
fig = Figure()
    ax11 = Axis(
        fig[1,1],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = L"c^{(1, 1)}"
        )
    ax12 = Axis(
        fig[1,2],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = L"c^{(1, 2)} = c^{(2, 1)}"
        )
    ax21 = Axis(
        fig[2,1],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = L"c^{(1, 2)} = c^{(2, 1)}"
        )

    ax22 = Axis(
        fig[2,2],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = L"c^{(2, 2)}"
        )
    boxplot!(ax11, ones(length(vec_cov11)),vec_cov11)
    boxplot!(ax12, ones(length(vec_cov12)),vec_cov12)
    boxplot!(ax21, ones(length(vec_cov21)),vec_cov21)
    boxplot!(ax22, ones(length(vec_cov22)),vec_cov22)
    linkyaxes!(ax11, ax12, ax21, ax22)
    fig
    save("/Users/berger/Documents/latex/pics//multidim_simu1_cov_box.png", fig)

# Efficiencies boxplots
fig = Figure()
    ax1 = Axis(
        fig[1,1],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = L"\rho_P^{(1)}"
        )
    ax2 = Axis(
        fig[1,2],
        xticklabelsvisible = false,
        xticksvisible = false,
        title = L"\rho_P^{(2)}"
        )
    boxplot!(ax1, ones(length(vec_eff1P)), vec_eff1P)
    boxplot!(ax2, ones(length(vec_eff2P)), vec_eff2P)
    linkyaxes!(ax1, ax2)
    fig
    save("/Users/berger/Documents/latex/pics//multidim_simu1_rho_box.png", fig)