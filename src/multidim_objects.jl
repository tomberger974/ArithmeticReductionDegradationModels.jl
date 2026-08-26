"""
    MaintenanceModels

Abstract type for Maintenance Models
"""
abstract type MaintenanceModels end


"""
    ARD1 submodel of the ARD models
"""
struct ARD1 <: MaintenanceModels end

"""
    ARDinf submodel of the ARD models
"""
struct ARDinf <: MaintenanceModels end


"""
    Efficiency

Efficiency object containing an imperfect maintenance model and a an efficiency value.

Need to be linked with a tuple (indicator, maintenance type) to have meaning
"""
mutable struct Efficiency
    value::Float64
    model::MaintenanceModels
end


"""
    MvWienerAR <: DPIM

Abstract type for Degradation Process with Imperfect Maintenance.
"""
abstract type DPIM end


"""
    MvWienerAR <: DPIM

Alias for the MvWienerAR type.
"""
const Indicator = Union{Symbol, String}
const EfficiencyKey = Tuple{Indicator, Indicator}
const Efficiencies = Dict{EfficiencyKey, Efficiency}

"""
    MvWienerAR <: DPIM

Alias for the MvWienerAR type.
"""
const CovarianceKey = Tuple{Indicator, Indicator}
const Covariances = Dict{CovarianceKey, Float64}


"""
    MvWienerAR <: DPIM

Multivariate (multiple indicators) Wiener degradation process with imperfect maintenance following a model of type AR.

Defined by a dictionary of drifts linking each indicator to its drift, a dictionary of volatilities linking each pair of indicators to their covariance, and a dictionary of efficiencies linking each pair (indicator, maintenance type) to an Efficiency object.
"""
mutable struct MvWienerAR <: DPIM
    drift::Dict{Indicator, Float64}
    covariances::Covariances
    efficiencies::Efficiencies

    function MvWienerAR(
            drift::AbstractDict,
            covariances::AbstractDict,
            efficiencies::AbstractDict)
        typed_drift = Dict{Indicator, Float64}()
        for (indicator, value) in drift
            indicator isa Indicator || throw(ArgumentError("drift keys must be Symbols or Strings"))
            typed_drift[indicator] = Float64(value)
        end

        typed_covariances = Dict{CovarianceKey, Float64}()
        for (key, value) in covariances
            key isa Tuple && length(key) == 2 && key[1] isa Indicator && key[2] isa Indicator || throw(ArgumentError("covariances keys must be pairs of Symbols or Strings"))
            typed_covariances[(key[1], key[2])] = Float64(value)
        end
        for ((ind1, ind2), value) in collect(typed_covariances)
            reverse_key = (ind2, ind1)
            if haskey(typed_covariances, reverse_key) && typed_covariances[reverse_key] != value
                throw(ArgumentError("covariances must be symmetric"))
            end
            typed_covariances[reverse_key] = value
        end

        typed_efficiencies = Dict{EfficiencyKey, Efficiency}()
        for (key, value) in efficiencies
            key isa Tuple && length(key) == 2 && key[1] isa Indicator && key[2] isa Indicator || throw(ArgumentError("efficiencies keys must be pairs of Symbols or Strings"))
            value isa Efficiency || throw(ArgumentError("efficiencies values must be Efficiency objects"))
            typed_efficiencies[(key[1], key[2])] = value
        end

        s = Set(keys(typed_drift))

        if Set([key[1] for key in keys(typed_covariances)]) != s || Set([key[2] for key in keys(typed_covariances)]) != s
            throw(ArgumentError("drift and covariances must share the same indicators"))
        end
        
        if Set([key[1] for key in keys(typed_efficiencies)]) != s
            throw(ArgumentError("drift and efficiencies must share the same indicators"))
        end

        new(typed_drift, typed_covariances, typed_efficiencies)
    end
end

MvWienerAR(; 
    drift::AbstractDict=Dict(:ind1 => 0.0),
    covariances::AbstractDict=Dict((:ind1, :ind1) => 1.0),
    efficiencies::AbstractDict=Dict((:ind1, :M) => Efficiency(.5, ARDinf()))) = MvWienerAR(drift, covariances, efficiencies)

MvWienerAR(indicators::Union{Vector{Symbol}, Vector{String}}, maintenance_types::Union{Vector{Symbol}, Vector{String}}) = MvWienerAR(Dict(ind => 0. for ind in indicators), Dict((ind1, ind2) => (ind1 == ind2 ? 1. : 0.) for ind1 in indicators, ind2 in indicators), Dict((ind, mt) => Efficiency(0., ARDinf()) for ind in indicators, mt in maintenance_types))


function Base.show(io::IO, mvw::MvWienerAR)
    print(io, "MvWienerAR(\n")

    print(io, "  drift:\n")
    for indicator in sort!(collect(keys(mvw.drift)), by=string)
        print(io, "    ", repr(indicator), " => ", mvw.drift[indicator], "\n")
    end

    print(io, "  covariances:\n")
    covariance_keys = sort!(collect(keys(mvw.covariances)), by=key -> (string(key[1]), string(key[2])))
    for key in covariance_keys
        print(io, "    ", repr(key), " => ", mvw.covariances[key], "\n")
    end

    print(io, "  efficiencies:\n")
    efficiency_keys = sort!(collect(keys(mvw.efficiencies)), by=key -> (string(key[1]), string(key[2])))
    for key in efficiency_keys
        efficiency = mvw.efficiencies[key]
        model_name = nameof(typeof(efficiency.model))
        print(io, "    ", repr(key), " => Efficiency(value = ", efficiency.value, ", model = ", model_name, "())\n")
    end

    print(io, ")")
end