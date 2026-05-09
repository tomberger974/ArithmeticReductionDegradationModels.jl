mutable struct MWARD1
    drift::Vector{Float64}
    volatility::Matrix{Float64}
    efficiencies::Dict{Symbol, Float64}
end

MWARD1(; drift::Vector{Float64}=[1.], volatility::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 1, 1), efficiencies::Dict{Symbol, Float64}=Dict(:M => 0.)) = MWARD1(drift, volatility, efficiencies)
MWARD1(dim::Int64; efficiencies::Dict{Symbol, Float64}=Dict(:M => 0.)) = MWARD1(ones(dim), diagm(ones(dim)), efficiencies)

Base.show(io::IO, mward1::MWARD1) = print(io, "drift:", mward1.drift, ", volatility:", mward1.volatility, ", efficiencies:", mward1.efficiencies)


"""
    Re-Simulate the values of the degradation according to `model` and modify consequently the degradations `VALUE` of the data set `data`.
    Optional parameter `from` indicates that the simulation of degradation `VALUE` is only done after the row indice `from` of `data`.
"""
function rand!(model::MWARD1, data::DegradationData, from::Int64=0)
#    ARDinf = (typeof(model) == WienerARD∞) #ARDinf or ARD1 ?
    modelParams = [model.drift, model.efficiencies]
    ρ = collect(values(model.efficiencies)) #initialisation
    deg = data.degradations
    maint = data.maintenances
    maintTypes = collect(keys(model.efficiencies))
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
            ρu = ρ[which_maint_type(maint.TYPE[i+1], maintTypes)]
        else
            t = vcat(t_maint, deg.DATE[i_maint])
            ip = 0
        end
        if (length(t) > 1)
            dt = t[2:end] - t[1:(end-1)]
            rd = rand(MvNormal(zeros(dt)), length(dt))
            #res[i] = rd
            degval = deg0 .+ cumsum(rd .* modelParams[2] .* sqrt.(dt) .+ modelParams[1] .* dt)
            deg.VALUE[i_maint] = degval[1:(end-ip)]
            if ARDinf
                deg0 = (1 - ρu) * degval[end]
            else
                #ARD1
                deg0 = degval[end] - ρu * (degval[end] - deg0)
            end
        end
    end

    data.infos["μ"] = modelParams[1]
    data.infos["σ"] = modelParams[2]
    data.infos["ρ"] = ρ
    data.infos["ρTypes"] = maintTypes
    if ARDinf
        data.infos["simulationModel"] = "ARDinf-Wiener"
    else
        data.infos["simulationModel"] = "ARD1-Wiener"
    end

    return data
end

function rand!(model::MWARD1, data::DegradationData)
    μ = model.drift
    Σ = model.volatility
    ρ = model.efficiencies
    deg = data.degradations
    maint = data.maintenances

    sort(vcat([0.], deg.DATE, maint.DATE))
    Δt = diff(sort(vcat([0.], deg.DATE, maint.DATE)))
    ΔN = rand.(MvNormal.([μ * x for x in Δt], [Σ * x for x in Δt]))
    W = reduce(hcat, vcat([[0., 0.]], cumsum(ΔN)))

    for i in 0:nrow(maint)
        if i == 0
            τi =  0.
            τiplus1 = maint[i+1]
            ρu = ρ[maint.TYPE[i+1]]
        elseif i == nrow(maint)
            τi = maint[i]
            τiplus1 = deg.DATE[end]
            ρu = ρ[maint.TYPE[i+1]]
        else
            τi = maint[i]
            τiplus1 = maint[i + 1]
            ρu = 1.
        end
        i_maint = findall(deg.NB_MAINTENANCES .== i)
        
    end

    return nothing
end