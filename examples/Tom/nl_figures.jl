using Latexify
using CairoMakie
using DataFrames
using Distributions

using Pkg
Pkg.activate(".")

using Revise
using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

k = 500
T = 1000
time = convert(Vector{Float64}, range(0, T, k))


#degradationdata parameters
k = 3#nb of maintenances
nj = 40#nb of obseravtions between each maintenances
T = (k+1)*nj
τ = zeros(k)
for i in 1:k
    τ[i] = i*(nj+1)
end
τ_types = ["P" for i in 1:k]
T = τ[end] + nj
Δt = 1.

#create the degradationdata instance
degradationdata = ARD.CreateData(τ, Δt, T, τ_types)
filter!(row -> row.TYPE in ["Between"], degradationdata.degradations)

nlwienerard1 = ARD.NLWienerARD1(1., .01, .5, 1., Dict("P" => .5))

rand!(nlwienerard1, degradationdata)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, deg.DATE, Y)
vlines!(ax, maint.DATE, color=:red)
display(fig)


Δtα = [ordered_time[i]^α - ordered_time[i-1]^α for i in 2:length(ordered_time)]
Δtβ = [ordered_time[i]^β - ordered_time[i-1]^β for i in 2:length(ordered_time)]

#creation of a sample of the distribution
rand(Normal(0, 1), 2)
Normal.(μ .* Δtα, σ .* sqrt.(Δtβ))
rand.(Normal.(μ .* Δtα, σ .* sqrt.(Δtβ)))
X = rand.(Normal.(μ .* Δtα, σ .* sqrt.(Δtβ)))