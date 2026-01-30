using ArithmeticReductionDegradationModels

#ARD1 parameters
μ = 2.
σ = sqrt(5)
ρ = [.5]

#create the wienerARD1 instance
wienerARD1 = WienerARD1(["P"], μ, σ, ρ)
wienerARD∞ = WienerARD∞(["P"], μ, σ, ρ)

#degradationdata parameters
k = 3#nb of maintenances
nj = 5#nb of obseravtions between each maintenances
T = (k+1)*nj
τ = zeros(k)
for i in 1:k
    τ[i] = i*(nj+1)
end
τ_types = ["P" for i in 1:k]
T = τ[end] + nj
Δt = 1.

#create the degradationdata instance
degradationdata = DegradationData(τ, Δt, T, τ_types)
filter!(row -> row.TYPE in ["Between"], degradationdata.degradations)
rand!(wienerARD1, degradationdata, 0)
fit_mle(wienerARD1, degradationdata)
fit_mle(wienerARD∞, degradationdata)