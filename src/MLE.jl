import Distributions.loglikelihood
import Distributions.fit_mle
import StatsBase.confint

"""
    computeSums(model::WienerARD∞, datas::Vector{DegradationData}, θs::Vector{Float64}; degTypes::Vector{String}=Vector{String}())

    Compute the different sums necessay to compute likelihood and MLE estimates.
     - `model` precise the `WienerARD∞` model used
     - `datas` a vector of iid DegradationData containing degradation observations and maintenance actions
     - `θs` the parameters value considered for computation
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes
     - `withDerivatives` specified if the first order and second order partial derivatives of the sums have also to be computed.
"""
function computeSums(model::WienerARD∞, datas::Vector{DegradationData}, θs::Vector{Float64}; degTypes::Vector{String}=Vector{String}(), withDerivatives::Bool=false)
    μ = θs[1]
    σ = θs[2]
    ρs = θs[3:end]
    nθ = length(θs)
    nρ = nθ - 2
    maintTypes = model.maintenances.TYPE

    #initialisations
    S1 = 0.
    S2 = 0.
    S3 = 0.
    S4 = 0.
    if withDerivatives
        #dρ
        dS1 = fill(0., nρ)
        d2S1 = fill(0., nρ, nρ)
        #dμ dρ
        dS2 = fill(0., nθ-1)
        d2S2 = fill(0., nθ-1, nθ-1)
    end
    N = 0.

    for l in 1:length(datas)
        data = datas[l]
        #println(string("data=", l))
        if isempty(degTypes)
            degradations = data.degradations
        else
            degradations = data.degradations[data.degradations.TYPE .∈ Ref(degTypes),:]
        end
        τs = vcat(0, data.maintenances.DATE, Inf64)
        types = vcat("missing", data.maintenances.TYPE, maintTypes[1])
        N += sum(unique(degradations.DATE) .> 0. )

        #Re-initialisation of computation variables for the l-th trajectory
        y_jmoins1_n_jmoins1 = 0.
        Aj = 0.
        Bj = 0.
        prod_j = 1.
        if withDerivatives
            #dρ
            dAj = fill(0., nρ)
            dBj = fill(0., nρ)
            dprod_j = fill(0., nρ)
            #dρ^2
            d2Aj = fill(0., nρ, nρ)
            d2Bj = fill(0., nρ, nρ)
            d2prod_j = fill(0., nρ, nρ)
        end

        for (j, τ) in enumerate(τs[1:(end-1)])
            degj = degradations[degradations.NB_MAINTENANCES .== j-1, :]
            #println(string("After τ_j=", τ))
            if isempty(degj.DATE)
                #no observation in the current maintenance interval
                U = which_maint_type(types[j+1], maintTypes)
                ρ = ρs[U]
                if withDerivatives
                    #dρ^2
                    d2Aj *= (1. - ρ)
                    d2Aj[:, U] -=  dAj
                    d2Aj[U, :] -=  dAj
                    d2Bj *= (1. - ρ)^2
                    d2Bj[:, U] -= 2 .* (1. - ρ) .* dBj
                    d2Bj[U, :] -= 2 .* (1. - ρ) .* dBj
                    d2Bj[U, U] += 2 .* ((τs[j+1] - τ) .+ Bj)
                    d2prod_j *= (1. - ρ) 
                    d2prod_j[:, U] -= dprod_j
                    d2prod_j[U, :] -= dprod_j
                    #dρ
                    dAj *= (1. - ρ)
                    dAj[U] -= (τs[j+1] - τ + Aj)
                    dBj *= (1. - ρ)^2 
                    dBj[U] -= 2 * (1. - ρ) * (τs[j+1] - τ + Bj)
                    dprod_j *= (1. - ρ) 
                    dprod_j[U] -= prod_j
                    #Be carrefull: computation of the derivative must be done before computation of the corresponding values since they are recursive
                end
                Aj = (1. - ρ) * (τs[j+1] - τ + Aj)
                Bj = (1. - ρ)^2 * (τs[j+1] - τ + Bj)
                prod_j *= (1. - ρ)
            else
                #observations in the current maintenance interval
                t_j_1, i = findmin(degj.DATE)
                y_j_1 = degj.VALUE[i]
                t_j_n_j, i = findmax(degj.DATE)
                y_j_n_j = degj.VALUE[i]
                Δt_j_1 = t_j_1 - τ
                Aj += Δt_j_1
                Bj += Δt_j_1
                #dρ and dρ^2 -> RAS: contributions do not depend on ρ
                #println(string("Aj=", Aj, " Bj=",Bj))
                if length(degj.DATE) > 1
                    if !issorted(degj.DATE)
                        ideg = sortperm(degj.DATE)
                        degj.DATE = degj.DATE[ideg]
                        degj.VALUE = degj.VALUE[ideg]
                    end
                    for ideg in 1:(length(degj.DATE)-1)
                        Δt_ji = degj.DATE[ideg+1] - degj.DATE[ideg]
                        S1 += log(Δt_ji)
                        S2 += (degj.VALUE[ideg+1] - degj.VALUE[ideg] - μ * Δt_ji )^2 / Δt_ji
                        if withDerivatives
                            #dμ
                            dS2[1] -= 2 *(degj.VALUE[ideg+1] - degj.VALUE[ideg] - μ * Δt_ji )
                            #dμ^2
                            d2S2[1,1] += 2 * Δt_ji
                            #dρ and dμ.dρ ->  RAS: contributions do not depend on ρ
                        end
                        #println(string("S2 +=", (degj.VALUE[ideg+1] - degj.VALUE[ideg] - μ * Δt_ji )^2 / Δt_ji, " (Between)"))
                    end
                end
                if Aj != 0.
                    S1 += log(Bj)
                    Cj = y_j_1 - prod_j * y_jmoins1_n_jmoins1 - μ * Aj #depend of μ
                    S2 += Cj^2 / Bj
                    S3 += y_j_n_j - y_j_1 + Aj / Bj * (y_j_1 - prod_j * y_jmoins1_n_jmoins1) 
                    S4 += t_j_n_j - t_j_1 + Aj^2 / Bj 
                    if withDerivatives
                        #dμ
                        dS2[1] -= 2 * Aj * Cj / Bj
                        #dρ
                        dS1 += dBj ./ Bj
                        dCj = dprod_j .* (- y_jmoins1_n_jmoins1) .- μ .* dAj #depend of μ
                        dS2[2:end] += dCj .* (2 * Cj / Bj) .- dBj .* (Cj^2 / Bj^2)
                        #dμ^2
                        d2S2[1,1] += 2 * Aj^2 / Bj
                        #dμ dρ
                        res = - dAj .* (2 * Cj / Bj) .+ dCj .* (-2 * Aj / Bj) .- dBj .* (-2 * Aj * Cj / Bj^2)
                        d2S2[2:end, 1] += res
                        d2S2[1, 2:end] += res
                        #dρ^2
                        d2S1 += d2Bj ./ Bj .- (dBj * transpose(dBj))./ Bj^2
                        d2Cj = d2prod_j .* (- y_jmoins1_n_jmoins1) .- μ .* d2Aj
                        d2S2[2:end, 2:end] += d2Cj .* (2 * Cj / Bj) .- d2Bj .* (Cj^2 / Bj^2) .+ 
                            (dCj * transpose(dCj .* (2 / Bj) .- dBj .* (2 * Cj / Bj^2) )) .- 
                            (dBj * transpose(dCj .* (2 * Cj / Bj^2) .- dBj .* (2 * Cj^2 / Bj^3)))
                    end
                    #println(string("S2 +=", (y_j_1 - (1. - ρ)^diff_ν_j * y_jmoins1_n_jmoins1 - μ * Aj)^2 / Bj, " (Befor)"))
                else
                    S3 += y_j_n_j - y_j_1 
                    S4 += t_j_n_j - t_j_1
                end 
                y_jmoins1_n_jmoins1 =  y_j_n_j
                U = which_maint_type(types[j+1], maintTypes)
                ρ = ρs[U]
                Aj = (1. - ρ) * (τs[j+1] - t_j_n_j)
                Bj = (1. - ρ)^2 * (τs[j+1] - t_j_n_j)
                prod_j = (1. - ρ)
                if withDerivatives
                    #dρ
                    dAj = fill(0., nρ)
                    dAj[U] = -(τs[j+1] - t_j_n_j)
                    dBj = fill(0., nρ)
                    dBj[U] = -2 * (1. - ρ) * (τs[j+1] - t_j_n_j)
                    dprod_j = fill(0., nρ)
                    dprod_j[U] = -1
                    #dρ^2
                    d2Aj = fill(0., nρ, nρ)
                    d2Bj = fill(0., nρ, nρ)
                    d2Bj[U,U] = 2 * (τs[j+1] - t_j_n_j)
                    d2prod_j = fill(0., nρ, nρ)
                end
            end
        end
    end

    if withDerivatives
        return (N, S1, S2, S3, S4, dS1, dS2, d2S1, d2S2) 
    else
        return (N, S1, S2, S3, S4)
    end
end

"""
    loglikelihood(model::WienerARD∞, datas::Vector{DegradationData}, θs::Vector{Float64}; degTypes::Vector{String}=Vector{String}(), withDerivatives::Bool=false)

    Compute the loglikelihood
     - `model` precise the `WienerARD∞` model used
     - `datas` a vector of iid DegradationData containing degradation observations and maintenance actions
     - `θs` the parameters value considered for computation
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes
     - `withDerivatives` specified if the first order and second order partial derivatives of the sums have also to be computed.
"""
function loglikelihood(model::WienerARD∞, datas::Vector{DegradationData}, θs::Vector{Float64}; degTypes::Vector{String}=Vector{String}(), withDerivatives::Bool=false)
    Ss = computeSums(model, datas, θs, degTypes=degTypes, withDerivatives=withDerivatives)
    σ = θs[2]
    nθ = length(θs)
    S1 = Ss[2]
    S2 = Ss[3]
    if withDerivatives
        #dρ
        dS1 = Ss[6]
        d2S1 = Ss[8]
        #dμ dρ
        dS2 = Ss[7]
        d2S2 = Ss[9]
    end
    N = Ss[1]
    
    lnL = - N/2 * log(2 * π ) - N * log(σ) - S1 / 2 - S2 / (2 * σ^2)
    if withDerivatives
        dlnL = Vector{Float64}(undef, nθ)
        d2lnL = Matrix{Float64}(undef, nθ, nθ)
        #dμ
        dlnL[1] = -dS2[1] / (2 * σ^2)
        #dσ
        dlnL[2] = - N / σ  + S2 / σ^3
        #dρ
        dlnL[3:end] .= -dS1 ./2 .- dS2[2:end] ./ (2 * σ^2)
        #dμ^2
        d2lnL[1,1] = -d2S2[1,1] / (2 * σ^2)
        #dμ dσ
        d2lnL[1,2] =  dS2[1] / σ^3
        d2lnL[2,1] =  d2lnL[1,2]
        #dσ^2
        d2lnL[2,2] =  N / σ^2  - 3 * S2 / σ^4
        #dμ dρ
        d2lnL[1,3:end] .= - d2S2[1,2:end] ./ (2 * σ^2)
        d2lnL[3:end,1] .= d2lnL[1,3:end]
        #dσ dρ
        d2lnL[2,3:end] .= dS2[2:end] ./  σ^3
        d2lnL[3:end,2] .= d2lnL[2,3:end]
        #dρ^2
        d2lnL[3:end,3:end] .= -d2S1 ./2 .- d2S2[2:end,2:end] ./ (2 * σ^2)
        return (lnL , dlnL, d2lnL)
    else
        return lnL
    end
end

"""
    profileloglikelihood(model::WienerARD∞, datas::Vector{DegradationData}, θs::Vector{Float64}; degTypes::Vector{String}=Vector{String}())

    Compute the profile loglikelihood, i.e. for a given value of the maintenance efficiency parameter, the value of the log-likelihood maximizes versus μ and σ
     - `model` precise the `WienerARD∞` model used
     - `datas` a vector of iid DegradationData containing degradation observations and maintenance actions
     - `ρs` the maintenance efficiency parameters value considered for computation
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes
"""
function profileloglikelihood(model::WienerARD∞, datas::Vector{DegradationData}, ρs::Vector{Float64}; degTypes::Vector{String}=Vector{String}())
    Ss = computeSums(model, datas, vcat(0, 1, ρs), degTypes=degTypes, withDerivatives=false)
    μ = Ss[4]/Ss[5]
    Ss = computeSums(model, datas, vcat(μ, 1, ρs), degTypes=degTypes, withDerivatives=false)
    N = Ss[1]
    S1 = Ss[2]
    S2 = Ss[3]
    σ = sqrt(S2 / N)
    #println(string("μ=",μ,"; σ=",σ,"; S1=",S1,"; S2=",S2,"; N=", N))
    #println(string("T1=", N * (log(2 * π)/2 + log(σ)), "; T2=", S1 / 2 ,"; T4=", N/2))
    return [ -( N * (log(2 * π)/2 + log(σ)) + S1 / 2 + N/2 ), μ , σ]    
end

# Internal function which verify if degradation is observed befor and after the same maintenance operation
# In fact, if this is the case, the exact value of ρ can be deduced from these observations
# In this case the function also verify that all the corresponding observations lead to the same value of ρ
# for the different considered maintenance types
# if this is not the case, there is a warning and ρ is computed according to the first observation corresponding to this maintenance type
# The function retruns the values of ρ if they can be such directly computed, 
# otherwise the function return missing for the maintenance types for which it can not be computed
function BeforeAndAfterMaintenance(model::WienerARD, datas::Vector{DegradationData}; degTypes::Vector{String}=Vector{String}())
    maintTypes = model.maintenances.TYPE
    if isempty(maintTypes)
        n = 1
    else
        n = length(maintTypes)
    end
    answear = Vector{Union{Missing, Float64}}(missing, n)
    first_l_answear = Vector{Union{Missing, Float64}}(missing, n)
    first_τ_answear = Vector{Union{Missing, Float64}}(missing, n)
    for l in 1:length(datas)
        data = datas[l]
        if isempty(degTypes)
            degradations = data.degradations
        else
            degradations = data.degradations[data.degradations.TYPE .∈ Ref(degTypes),:]
        end
        for (j,τ) in enumerate(data.maintenances.DATE)
            degτ = degradations[degradations.DATE .== τ, :]
            if nrow(degτ) > 1
                i = which_maint_type(data.maintenances.TYPE[j], maintTypes)
                if ismissing(answear[i])
                    answear[i] = 1 - degτ.VALEUR[degτ.NB_MAINTENANCES .== j][1] /  degτ.VALUE[degτ.NB_MAINTENANCES .== j-1][1]
                    first_l_answear[i] = l
                    first_τ_answear[i] = τ
                elseif answear[i] != 1 - degτ.VALEUR[degτ.NB_MAINTENANCES .== j][1] /  degτ.VALUE[degτ.NB_MAINTENANCES .== j-1][1]
                    @warn string("The model is not coherent with the observations of degradation befor and after maintenance between systems ",first_l_answear, " and ", l, " and between maintenance times ", first_τ_answear, " and ",τ)
                end
            end
        end
    end
    return answear
end

"""
    function fit_mle(model::WienerARD∞, datas::Vector{DagradationData}; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)

    Compute the maximum likelihood estimators
     - `model` precise the `WienerARD∞` model used
     - `datas` a vector of iid DegradationData containing degradation observations and maintenance actions
     - `ρ0s` the initialisation value of maintenance efficiency parameters (not the underlying degradation paramters) for profile likelihood maximization
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes
"""
function fit_mle(model::WienerARD∞, datas::Vector{DegradationData}; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)
    ρs = BeforeAndAfterMaintenance(model, datas, degTypes=degTypes)
    ρsMissing = ismissing.(ρs)
    ρsobserved = fill(false, length(ρs))#which maintenance type is effectively observed on datas
    for j in 1:length(datas)
        datasimaintType = unique(datas[j].maintenances.TYPE)
        ρsobserved = ρsobserved .|| [model.maintenances.TYPE[i] ∈ datasimaintType for i in 1:length(model.maintenances.TYPE)]
    end
    ρs[ρsMissing] .= 0
    ρs = convert(Vector{Float64}, ρs)
    if sum(ρsMissing) > 0
        res_optim = optimize(x-> (ρs[ρsMissing .&& ρsobserved] .= x; -(profileloglikelihood(model, datas, ρs; degTypes=degTypes))[1]), ρ0s[ρsMissing .&& ρsobserved])
        if print_summary
            println(summary(res_optim))
        end
        ρs[ρsMissing .&& ρsobserved] .= Optim.minimizer(res_optim)
    end
    res = profileloglikelihood(model, datas, ρs, degTypes=degTypes)
    if sum(.!ρsobserved) > 0
        ρs = convert(Vector{Union{Float64, Missing}}, ρs)
        ρs[.!ρsobserved] .= missing
    end
    estim = DataFrame(transpose(vcat(res[2:3], ρs)), :auto)
    rename!(estim, vcat(["μ", "σ"], string.("ρ_", model.maintenances.TYPE)))
    return estim
end

"""
    fit_mle(model::WienerARD∞, data::DegradationData; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)

    Compute the maximum likelihood estimators
     - `model` precise the `WienerARD∞` model used
     - `data` a DegradationData containing degradation observations and maintenance actions
     - `ρ0s` the initialisation value of maintenance efficiency parameters (not the underlying degradation paramters) for profile likelihood maximization
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes

"""
function fit_mle(model::WienerARD∞, data::DegradationData, ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)); degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)
    return fit_mle(model, [data], ρ0s=ρ0s, degTypes=degTypes, print_summary=print_summary)
end


"""
    function fit_mle!(model::WienerARD∞, datas::Vector{DagradationData}; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)

    Update the `model` parameters value with the corresponding computed maximum likelihood estimators
     - `model` precise the `WienerARD∞` model used
     - `datas` a vector of iid DegradationData containing degradation observations and maintenance actions
     - `ρ0s` the initialisation value of maintenance efficiency parameters (not the underlying degradation paramters) for profile likelihood maximization
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes
"""
function fit_mle!(model::WienerARD∞, datas::Vector{DegradationData}; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)
    estim = collect(fit_mle(model, datas, ρ0s=ρ0s, degTypes=degTypes, print_summary=print_summary)[1,:])
    model.underlyingDegradation = estim[1:2]
    for (i, ρ) in enumerate(estim[3:end])
        if ismissing(ρ)
            model.maintenances.VALUE[i] = ρ0s[i]
        else
            model.maintenances.VALUE[i] = ρ
        end
    end
    return model
end

"""
    fit_mle!(model::WienerARD∞, data::DegradationData; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)

    Update the `model` parameters value with the corresponding computed maximum likelihood estimators
     - `model` precise the `WienerARD∞` model used
     - `data` a DegradationData containing degradation observations and maintenance actions
     - `ρ0s` the initialisation value of maintenance efficiency parameters (not the underlying degradation paramters) for profile likelihood maximization
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes

"""
function fit_mle!(model::WienerARD∞, data::DegradationData, ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)); degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)
    return fit_mle!(model, [data], ρ0s=ρ0s, degTypes=degTypes, print_summary=print_summary)
end

"""
    confint(model::WienerARD∞, datas::Vector{DegradationData}, α::Float64=0.05; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)

    Compute the asymptotic confidence interval for parameters maximum likelihood estimation based on estimated Fisher information
     - `model` precise the `WienerARD∞` model used
     - `datas` a vector of iid DegradationData containing degradation observations and maintenance actions
     - `level` the confidence interval level 
     - `ρ0s` the initialisation value of maintenance efficiency parameters (not the underlying degradation paramters) for profile likelihood maximization
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes

"""
function confint(model::WienerARD∞, datas::Vector{DegradationData}, level::Float64=0.95; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)
    θest = fit_mle(model, datas, ρ0s=ρ0s, degTypes=degTypes, print_summary=print_summary)
    missingEst = ismissing.(collect(θest[1,:]))
    θest[1, missingEst] = ρ0s[missingEst[3:end]]
    α = 1. - level
    u = quantile(Normal(),1-α/2)
    In = loglikelihood(model, datas, collect(θest[1,:]), withDerivatives=true)[3]
    In = In[.!missingEst, .!missingEst]
    sqrt_inv_In = - Matrix(1.0I, size(In)[1], size(In)[2]) / In
    sqrt_inv_In = sqrt(sqrt_inv_In )
    if typeof(sqrt_inv_In) == Matrix{ComplexF64}
        @warn("Conditionning problem with the loglikelihood hessian computed at the estimated parameter values.")
        return DataFrame(Matrix{Union{Float64, Missing}}(missing, 2, ncol(θest)), :auto)
    else
        if sum(missingEst)>0
            intConfAss = DataFrame(Matrix{Union{Float64, Missing}}(missing, 2, ncol(θest)), :auto)
        else
            intConfAss = DataFrame(Matrix{Float64}(undef, 2, ncol(θest)), :auto)
        end
        rename!(intConfAss, names(θest))
        intConfAss[1, .!missingEst] = collect(θest[1,.!missingEst]) .- u .* [sqrt_inv_In[i,i] for i in 1:size(sqrt_inv_In)[1]]
        intConfAss[2, .!missingEst] = collect(θest[1,.!missingEst]) .+ u .* [sqrt_inv_In[i,i] for i in 1:size(sqrt_inv_In)[1]]
        return intConfAss
    end
end

"""
    confint(model::WienerARD∞, data::DegradationData, level::Float64=0.95; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)

    Compute the asymptotic confidence interval for parameters maximum likelihood estimation based on estimated Fisher information
     - `model` precise the `WienerARD∞` model used
     - `datas` a DegradationData containing degradation observations and maintenance actions
     - `level` the confidence interval level 
     - `ρ0s` the initialisation value of maintenance efficiency parameters (not the underlying degradation paramters) for profile likelihood maximization
     - `degTypes` can possibly specified that only some degradation types are considered, i.e. the degradation observations forwhich the value in `TYPE` column is in degTypes

"""
function confint(model::WienerARD∞, data::DegradationData, level::Float64=0.95; ρ0s::Vector{Float64}=fill(0.5, nrow(model.maintenances)), degTypes::Vector{String}=Vector{String}(), print_summary::Bool=true)
    return confint(model, [data], level, ρ0s=ρ0s, degTypes=degTypes, print_summary=print_summary)
end