using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

using DataFrames
using Distributions
using LinearAlgebra

ARD.MWARD1(drift=[2.])
mward1 = ARD.MWARD1(2; efficiencies=Dict(:P => 1., :C => 0.))
μ = mward1.drift
Σ = mward1.volatility
ρ = mward1.efficiencies

collect(keys(mward1.efficiencies))

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
maint = degradationdata.maintenances
deg = degradationdata.degradations

t = vcat([0.], maint.DATE, deg.DATE)
p = sortperm(t)
Δt = diff(sort(vcat([0.], deg.DATE, maint.DATE)))
Δμ = [μ * x for x in Δt]
ΔΣ = [Σ * x for x in Δt]
ΔN = rand.(MvNormal.(Δμ, ΔΣ))
W = reduce(hcat, vcat([[0., 0.]], cumsum(ΔN)))

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, sort(vcat([0.], deg.DATE, maint.DATE)), W[1, :])
display(fig)

vcat(0., maint, )



rand(MvNormal(zeros(2), diagm(ones(2))))
rand!(mward1, degradationdata, 0)