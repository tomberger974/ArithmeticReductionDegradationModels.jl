
import Random.rand!, Base.rand

"""
    rand!(model::WienerARD, data::DegradationData, from::Int64=0)

    Re-Simulate the values of the degradation according to `model` and modify consequently the degradations `VALUE` of the data set `data`.
    Optional parameter `from` indicates that the simulation of degradation `VALUE` is only done after the row indice `from` of `data`.

"""
function rand!(model::WienerARD, data::DegradationData, from::Int64=0, drift::Function=identity, volatility::Function=identity)
    try
        # Attempt to call `drift` with Float64 input
        result = drift(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Drift function must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `drift` failed, clarify the error
        throw(ArgumentError("Drift function must accept Float64 input"))
    end

    try
        # Attempt to call `volatility` with Float64 input
        result = volatility(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Volatility function must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `volatility` failed, clarify the error
        throw(ArgumentError("Volatility function must accept Float64 input"))
    end

    ARDinf = (typeof(model) == WienerARD∞) #ARDinf or ARD1 ?
    modelParams = collect(params(model)[1,:])
    ρ = modelParams[3] #initialisation
    deg = data.degradations
    maint = data.maintenances
    maintTypes = model.maintenances.TYPE
    if from == 0
        t0 = 0.
        deg0 = 0.
        i0 = 0
    else
        t0 = deg.DATE[min(from, nrow(data.degradations))]
        deg0 = deg.VALUE[min(from, nrow(data.degradations))]
        i0 = deg.NB_MAINTENANCES[min(from, nrow(data.degradations))]
    end
    for i in i0:nrow(maint)
        t_maint = max( (vcat(0., maint.DATE))[i+1] , t0 )
        i_maint = findall(deg.NB_MAINTENANCES .== i)
        i_maint = i_maint[i_maint .> from]
        if i< length(maint.DATE)
            t = vcat(t_maint, deg.DATE[i_maint],maint.DATE[i+1])
            ip = 1
            ρ = modelParams[2 + which_maint_type(maint.TYPE[i+1], maintTypes)]
        else
            t = vcat(t_maint, deg.DATE[i_maint])
            ip = 0
        end
        if (length(t) > 1)
            dtdrift = [drift(x) for x in t[2:end]] - [drift(x) for x in t[1:(end-1)]]
            dtvolatility = [volatility(x) for x in t[2:end]] - [volatility(x) for x in t[1:(end-1)]]
            rd = rand(Normal(), length(dtdrift))
            #res[i] = rd
                degval = deg0 .+ cumsum(rd .* modelParams[2] .* sqrt.(dtvolatility) .+ modelParams[1] .* dtdrift)
            deg.VALUE[i_maint] = degval[1:(end-ip)]
            if ARDinf
                deg0 = (1 - ρ) * degval[end]
            else
                #ARD1
                deg0 = degval[end] - ρ * (degval[end] - deg0)
            end
        end
    end

    data.infos["μ"] = modelParams[1]
    data.infos["σ"] = modelParams[2]
    data.infos["ρ"] = modelParams[3:end]
    data.infos["ρTypes"] = maintTypes
    if ARDinf
        data.infos["simulationModel"] = "ARDinf-Wiener"
    else
        data.infos["simulationModel"] = "ARD1-Wiener"
    end

    return data
end


# Create a new SystData object and simulate its degradation values
# for special case of periodic maintenances and observations
# Δτ: time between maintenances
# nbτ: number of maintenances
# Δt: time between degradation observations
# params = (μ, σ, ρ)
# Degradation observations have a TYPE specifying if there are done "After", "Before" or "Between" maintenances

"""
    Base.rand(model::WienerARD, Δτ::Float64, nbτ::Int64, Δt::Float64)

    Create a degradation data set `DegradationData` and simulate its degradation values for special case of
        - `nbτ` periodic maintenances every `Δτ` time units,
        -  periodic observations every `Δτ` up to `Δτ * (nbτ+1)`.
    Degradation observations have a TYPE specifying if there are done "After", "Before" or "Between" maintenances

"""
function rand(model::WienerARD, Δτ::Float64, nbτ::Int64, Δt::Float64)
    return rand(model, convert(Vector{Float64}, Δτ .* (1 : nbτ)), Δt, Δτ * (nbτ+1))
end

"""
    Base.rand(model::WienerARD, maintenancesDates::Vector{Float64}, Δt::Float64, t_max::Float64=maintenancesDates[end], maintenancesTypes::Vector{Float64}=fill("",length(maintenancesDates)))

    Create a degradation data set `DegradationData` and simulate its degradation values for special case of
        - maintenances at times `maintenancesDates` with types that can be specifed by `maintenancesTypes`
        - periodic degradation observations every `Δt` units of time up to time specified optionally by `t_max` or corresponding to the last maintenance time
    Degradation observations have a TYPE specifying if there are done "After", "Before" or "Between" maintenances

"""
function rand(model::WienerARD, maintenanceDates::Vector{Float64}, Δt::Float64, t_max::Float64=maintenanceDates[end], maintenanceTypes::Vector{String}=fill("",length(maintenanceDates)))

    data = CreateData(maintenanceDates, Δt, t_max, maintenanceTypes) # function CreateData in Tools.jl
    rand!(model, data)
    return data
end

"""
    rand!(model::WienerARD, datas::Vecor{DegradationData}, froms::Vecor{Int64}=fill(0, length(datas)))

    Re-Simulate the values of the degradation according to `model` and modify consequently the degradations `VALUE` of the data set `datas`.
    Optional parameter `froms` indicates that the simulation of degradation `VALUE` is only done after the row indice `from[i]` for `datas[i]`.
    Be carreful `froms` and `datas` must have the same length
"""
function rand!(model::WienerARD, datas::Vector{DegradationData}, froms::Vector{Int64}=fill(0, length(datas)))
    if length(datas) != length(froms)
        error("Arguments datas and froms must have the same length") 
    end
    for i in 1:length(datas)
        rand!(model, datas[i], froms[i])
    end
    return datas
end