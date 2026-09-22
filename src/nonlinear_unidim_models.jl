abstract type NLWienerARD end

mutable struct NLWienerARD1 <: NLWienerARD
    drift::Float64
    dispertion::Float64
    power_drift::Float64
    power_dispersion::Float64
    efficiencies::Dict{Union{Symbol, String}, Float64}
end

NLWienerARD1(; μ::Float64=1., σ::Float64=1., α::Float64=1., β::Float64=1., ρ = Dict(:M => 0.)) = NLWienerARD1(μ, σ, α, β, ρ)

mutable struct NLWienerARD∞ <: NLWienerARD
    drift::Float64
    dispertion::Float64
    power_drift::Float64
    power_dispersion::Float64
    efficiencies::Dict{Union{Symbol, String}, Float64}
end

NLWienerARD∞(; μ::Float64=1., σ::Float64=1., α::Float64=1., β::Float64=1., ρ = Dict(:M => 0.)) = NLWienerARD∞(μ, σ, α, β, ρ)