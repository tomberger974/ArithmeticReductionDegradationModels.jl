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
K = 3 #number of maintenances
N_i = 5 #number of observations between each maintenance
Δt = 1. #time increment between 2 observations
τ = convert(Vector{Float64}, [j*Δt*N_i for j in 1:K]) #maintenance dates
τ_types = ["P" for _ in 1:K] #maintenance types
T = convert(Float64, τ[end] + Δt*N_i) #maximal time of observation


#create the degradationdata instance
degradations = DataFrame(DATE = 1:Δt:T, VALUE1 = Vector{Float64}(undef, Int64(T/Δt)), VALUE2 = Vector{Float64}(undef, Int64(T/Δt)), NB_MAINTENANCES = count_NB_MAINTENANCES(τ, Vector{Float64}(1:Δt:T)))
maintenances = DataFrame(DATE = τ, TYPE = τ_types)
degradationdata = ARD.MvDegradationData(maintenances, degradations)
filter!(row -> row.DATE ∉ degradationdata.maintenances.DATE, degradationdata.degradations)
ARD.rand!(mward1, degradationdata)


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

