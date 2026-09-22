import Distributions.params
import Distributions.mean
import Distributions.var
import Distributions.std
import Distributions.quantile

abstract type WienerARD end

"""
    Define a uni-dimensional ARD∞ (Arithmtic Reduction of the Degradation) model with a underlying degradation following a Wiener process. 
    ARD∞ assumption is that maintenances reduce the degradation proportionnally to its value just before maintenance.
"""
mutable struct WienerARD∞ <: WienerARD
    underlyingDegradation::Vector{Float64}
    maintenances::DataFrame

    """
    WienerARD∞(maintTypes::Vector{String}, μ::Float64 = 0., σ::Float64 = 1., ρs::Vector{Float64} = (length(maintTypes) == 0 ? convert(Vector{Float64},[]) : fill(0.5, length(maintTypes))))

    Constructor function for ARD∞ model type, where
        `μ` : the linear slop parameter of the underlying degradation
        `σ` : the dispertion parameter of the underlying degradation
        `ρs` : the possibly multidimentional maintenance efficiency parameter
        `maintTypes` : the corresponding mainteance types identifier
    """
    function WienerARD∞(maintTypes::Vector{String}, μ::Float64 = 0., σ::Float64 = 1., ρs::Vector{Float64} = (length(maintTypes) == 0 ? convert(Vector{Float64},[]) : fill(0.5, length(maintTypes))))
        if length(ρs) != length(maintTypes)
            error("The vectors of maintence parameters value and maintenances type must have the same length.")
        end
        new([μ, σ], DataFrame(TYPE = maintTypes, VALUE = ρs))
    end

    """
    WienerARD∞(maintType::String = "", μ::Float64=0., σ::Float64=1., ρ::Float64=0.5)

    Constructor function for ARD∞ model type, where 
        `μ` : the linear slop parameter of the underlying degradation
        `σ` : the dispertion parameter of the underlying degradation
        `ρ` : a single maintenance efficiency parameter
        `maintType` : a single maintenance type identifier
    """
    function WienerARD∞(maintType::String = "", μ::Float64 = 0., σ::Float64 = 1., ρ::Float64 = 0.5)
        new([μ, σ], DataFrame(TYPE = [maintType], VALUE = [ρ]))
    end
end

"""
    Define a uni-dimensional ARD∞ (Arithmtic Reduction of the Degradation) model with a underlying degradation following a Wiener process. 
    ARD∞ assumption is that maintenances reduce the degradation proportionnally to its value just before maintenance.
"""
mutable struct WienerARD1 <: WienerARD
    underlyingDegradation::Vector{Float64}
    maintenances::DataFrame

    """
    WienerARD1(maintTypes::Vector{String}, μ::Float64 = 0., σ::Float64 = 1., ρ::Vector{Float64} = (length(maintTypes) == 0 ? convert(Vector{Float64},[]) : fill(0.5, length(maintTypes))))

    Constructor function for ARD1 model type, where
        `μ` : the linear slop parameter of the underlying degradation
        `σ` : the dispertion parameter of the underlying degradation
        `ρs` : the possibly multidimentional maintenance efficiency parameter
        `maintTypes` : the corresponding mainteance types identifier
    """
    function WienerARD1(maintTypes::Vector{String}, μ::Float64 = 0., σ::Float64 = 1., ρ::Vector{Float64} = (length(maintTypes) == 0 ? convert(Vector{Float64},[]) : fill(0.5, length(maintTypes))))
        if length(ρ) != length(maintTypes)
            error("The vectors of maintence parameters value and maintenances type must have the same length.")
        end
        new([μ, σ], DataFrame(TYPE = maintTypes, VALUE = ρ))
    end


    """
    WienerARD1(maintType::String = "", μ::Float64 = 0., σ::Float64 = 1., ρ::Float64 = 0.5)

    Constructor function for ARD1 model type, where 
        `μ` : the linear slop parameter of the underlying degradation
        `σ` : the dispertion parameter of the underlying degradation
        `ρ` : a single maintenance efficiency parameter
        `maintType` : a single maintenance type identifier
    """
    function WienerARD1(maintType::String = "", μ::Float64 = 0., σ::Float64 = 1., ρ::Float64 = 0.5)
        new([μ, σ], DataFrame(TYPE = [maintType], VALUE = [ρ]))
    end
end

"""
    Distributions.params(model::WienerARD)

Return a DataFrame specifyng the parameters of the ARD model.
"""
function params(model::WienerARD)
    p = DataFrame(transpose(vcat(model.underlyingDegradation, model.maintenances.VALUE)), :auto)
    nn = Vector{String}(undef, nrow(model.maintenances) + 2)
    nn[1:2] .= ["μ", "σ"]
    for i in 1:nrow(model.maintenances)
        if model.maintenances.TYPE[i] == ""
            nn[i + 2] = "ρ"
        else
            nn[i + 2] = string("ρ_", model.maintenances.TYPE[i])
        end
    end
    rename!(p, nn)
    return p
end


"""
    params!(model::WienerARD, θ::Vector{Float64})

Modify the parameters value of `model`.
"""
function params!(model::WienerARD, θ::Vector{Float64})
    if length(θ) != nrow(model.maintenances) + 2
        error("The dimension of the parameter vector is not coherent with the model.")
    end
    model.underlyingDegradation .= θ[1:2]
    model.maintenances.VALUE .= θ[3:end]
    return params(model)
end


"""
    params!(model::WienerARD, θ::Vector{Float64}, maintTypes::Vector{String})

Modify the parameters value and maintenance types of `model`.
"""
function params!(model::WienerARD, θ::Vector{Float64}, maintTypes::Vector{String})
    if length(θ) != nrow(model.maintenances) + 2
        error("The dimension of the parameter vector is not coherent with the model.")
    end
    if length(maintTypes) != nrow(model.maintenances)
        error("The dimension of the maintenance types vector is not coherent with the model.")
    end
    model.underlyingDegradation .= θ[1:2]
    model.maintenances.VALUE .= θ[3:end]
    model.maintenances.TYPE .= maintTypes
    return params(model)
end


"""
    mean_var(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))

    Compute the mean and variance of a WienerARD∞ `model` at the time at which degradation is observed in `data` and following also the maintenance times and types of `data``.
    Optional argument `from` and `to` corresponds to rows indexes in the degradation DataFrame `data`: mean and variance are computed conditionnaly to the degradation observation value at index `from`  and up to the degradation observation of index `to` 
"""
function mean_var(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))
    θ = collect(params(model)[1,:])
    maintTypes = model.maintenances.TYPE

    if from > to 
        error("Argument from must be greater than argument to.")
    end

    to = min(to, nrow(data.degradations))

    res = DataFrame(
        DATE = vcat(0, data.degradations.DATE[1:to]),
        NB_MAINTENANCES = vcat(0, data.degradations.NB_MAINTENANCES[1:to]),
        MEAN= zeros(Float64, to+1),
        VAR = zeros(Float64, to+1)
        )
    if from == 0
        nb_PM = 0
    else
        from = min(from, to+1)
        for i in 1:from
            res.MEAN[i+1] = data.degradations.VALUE[i]
        end      
        nb_PM = data.degradations.NB_MAINTENANCES[from]
    end
    for i in (from+1):to
        #Notice that : res.DATE[i+1] = data.degradations.DATE[i]
        #The interest of res.DATE is that it also integrate time 0

        if nb_PM < data.degradations.NB_MAINTENANCES[i]
            #at least 1 maintenance since last observation
            nb_PM += 1
            ρ = θ[2 + which_maint_type(data.maintenances.TYPE[nb_PM], maintTypes)]
            res.MEAN[i+1] = (1-ρ) * (res.MEAN[i] + θ[1] * (data.maintenances.DATE[nb_PM] - res.DATE[i]))
            res.VAR[i+1] = (1-ρ)^2 * (res.VAR[i] + θ[2]^2 * (data.maintenances.DATE[nb_PM] - res.DATE[i]))
            while nb_PM < data.degradations.NB_MAINTENANCES[i]
                # more than 1 maintenance since the last observation
                nb_PM += 1
                ρ = θ[2 + which_maint_type(data.maintenances.TYPE[nb_PM], maintTypes)]
                res.MEAN[i+1] = (1-ρ) * (res.MEAN[i+1] + θ[1] *(data.maintenances.DATE[nb_PM] - data.maintenances.DATE[nb_PM-1]))
                res.VAR[i+1] = (1-ρ)^2 * (res.VAR[i+1] + θ[2]^2 *(data.maintenances.DATE[nb_PM] - data.maintenances.DATE[nb_PM-1]))
            end
            res.MEAN[i+1] = res.MEAN[i+1]  + θ[1] * (res.DATE[i+1] - data.maintenances.DATE[nb_PM])
            res.VAR[i+1] = res.VAR[i+1] + θ[2]^2 * (res.DATE[i+1] - data.maintenances.DATE[nb_PM])
        else
            #no maintenance since last observation
            res.MEAN[i+1] = res.MEAN[i] +   θ[1] * (res.DATE[i+1] - res.DATE[i])
            res.VAR[i+1] = res.VAR[i] +  θ[2]^2 * (res.DATE[i+1] - res.DATE[i])
        end
    end
    return res[2:nrow(res),1:ncol(res)]
end

"""
    mean(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))

    Compute the mean of a WienerARD∞ `model` at the time at which degradation is observed in `data` and following also the maintenance times and types of `data``.
    Optional argument `from` and `to` corresponds to rows indexes in the degradation DataFrame `data`: mean is computed conditionnaly to the degradation observation value at index `from`  and up to the degradation observation of index `to` 
"""
function mean(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))
    return mean_var(model, data; from=from, to=to)[!,["DATE", "NB_MAINTENANCES", "MEAN"]]
end

"""
    var(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))

    Compute the variance of a WienerARD∞ `model` at the time at which degradation is observed in `data` and following also the maintenance times and types of `data``.
    Optional argument `from` and `to` corresponds to rows indexes in the degradation DataFrame `data`: variance is computed conditionnaly to the degradation observation value at index `from`  and up to the degradation observation of index `to` 
"""
function var(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))
    return mean_var(model, data; from=from, to=to)[!,["DATE", "NB_MAINTENANCES", "VAR"]]
end

"""
    sd(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))

    Compute the standard deviation of a WienerARD∞ `model` at the time at which degradation is observed in `data` and following also the maintenance times and types of `data``.
    Optional argument `from` and `to` corresponds to rows indexes in the degradation DataFrame `data`: standard deviation is computed conditionnaly to the degradation observation value at index `from`  and up to the degradation observation of index `to` 
"""
function std(model::WienerARD∞, data::DegradationData; from::Int64=0, to::Int64=typemax(Int64))
    res = mean_var(model, data; from=from, to=to)[!, ["DATE", "NB_MAINTENANCES", "VAR"]]
    res[!, "STD"] = sqrt.(res[!, "VAR"])
    return res[!, ["DATE", "STD"]]
end

"""
    quantile(model::WienerARD∞, data::DegradationData, probs::Vecor{Float64}=[0.5]; from::Int64=0, to::Int64=typemax(Int64))

    Compute the quantile of a WienerARD∞ `model` at the time at which degradation is observed in `data` and following also the maintenance times and types of `data``.
    Optional argument `from` and `to` corresponds to rows indexes in the degradation DataFrame `data`: standard deviation is computed conditionnaly to the degradation observation value at index `from`  and up to the degradation observation of index `to` 
"""
function quantile(model::WienerARD∞, data::DegradationData, probs::Vector{Float64}=[0.5]; from::Int64=0, to::Int64=typemax(Int64))
    res = mean_var(model, data; from=from, to=to)
    quant = Vector{Float64}(undef, nrow(res))
    for prob in probs
        for i in 1:nrow(res) 
            quant[i] = quantile(Normal(res[i,"MEAN"], sqrt(res[i,"VAR"])), prob)
        end
        res[!, string("Q_",prob)] .= quant

    end
    return res[!, vcat(["DATE", "NB_MAINTENANCES"], string.("Q_",probs))]
end