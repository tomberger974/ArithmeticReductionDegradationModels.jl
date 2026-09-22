using CairoMakie
using Distributions

#####GENERATEUR DE DONNEES ALEATOIRES NON LINEAIRES#####
function generateur(wienerARD1::WienerARD1, degradationdata::DegradationData, f::Function=identity, g::Function=identity)
    try
        # Attempt to call `f` with Float64 input
        result = f(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Drift function `f` must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `f` failed, clarify the error
        throw(ArgumentError("Drift function `f` must accept Float64 input"))
    end

    try
        # Attempt to call `f` with Float64 input
        result = g(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Dispersion function `g` must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `f` failed, clarify the error
        throw(ArgumentError("Dispersion function `g` must accept Float64 input"))
    end

    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    μ = wienerARD1.underlyingDegradation[1]
    σ = wienerARD1.underlyingDegradation[2]
    type_to_ρ = Dict(wienerARD1.maintenances.TYPE .=> wienerARD1.maintenances.VALUE)

    NBM = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
    ν = NBM .+ 1
    K = length(ν)

    #Y2 takes the successive simulated values of the degradation level over time
    #Y1 corresponds to Y(τ_j+) which is necessary to compute Y(τ_{j+1}+) = Y(τ_{j+1}-) - ρ[Y(τ_{j+1}-) - Y(τ_j+)] 
    Y1 = 0.
    Y2 = 0.

    new_deg = Vector{Float64}(undef, 1)#the vector which will contain the new values for degradationdataz.degradations.VALUE

    #The first increment is special 
    if ν[1] == 1
        Δft = f(deg[1, "DATE"])
        Δgt = g(deg[1, "DATE"])
        Y2 += rand(Normal(μ * Δft, σ * sqrt(Δgt)))
        new_deg[1] = Y2
    else
        Δfτ = f(maint[1, "DATE"])
        Δgτ = g(maint[1, "DATE"])
        Y2 += rand(Normal(μ * Δfτ, σ * sqrt(Δgτ)))
        Y1 = Y2 - type_to_ρ[maint[1, "TYPE"]] * (Y2 - Y1)

        if ν[1] != 2
            for j in 2:ν[1]-1
                Δfτ = f(maint[j, "DATE"]) - f(maint[j-1, "DATE"])
                Δgτ = g(maint[j, "DATE"]) - g(maint[j-1, "DATE"])
                Y2 = Y1 + rand(Normal(μ * Δfτ, σ * sqrt(Δgτ)))
                Y1 = Y2 - type_to_ρ[maint[j, "TYPE"]] * (Y2 - Y1)
            end
        end

        Δft = f(deg[1, "DATE"]) - f(maint[deg[1, "NB_MAINTENANCES"], "DATE"])
        Δgt = g(deg[1, "DATE"]) - g(maint[deg[1, "NB_MAINTENANCES"], "DATE"])
        Y2 = Y1 + rand(Normal(μ * Δft, σ * sqrt(Δgt)))
        new_deg[1] = Y2
    end

    #For each ν[i], we generate y_{ν_i, j} for j≥2 and then generate y_{ν_{i+1}, 1} (if ν_{i+1} exists)
    for i in 1:K
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i] - 1, deg)#keep only the observations between \nu[i]-1 and \nu[i]

        if nrow(ν_ith_deg) != 1
            Δft = diff([f(t) for t in ν_ith_deg.DATE])
            Δgt = diff([g(t) for t in ν_ith_deg.DATE])
            ν_ith_deg[2:end, "VALUE"] = new_deg[end] .+ cumsum(rand.(Normal.(μ .* Δft, σ .* sqrt.(Δgt))))
        end

        if i != lastindex(ν)
            u_ν_i = maint[ν[i], "TYPE"]
            ρ_u_ν_i = type_to_ρ[u_ν_i]
            
            idx = findfirst(deg.NB_MAINTENANCES .== ν[i + 1] - 1)
            ν_iplus1th_first_deg = deg[[idx], :]#find the first observation following all the observations stored in ν_ith_deg

            Δft = f(maint[ν[i], "DATE"]) - f(ν_ith_deg[end, "DATE"])
            Δgt = g(maint[ν[i], "DATE"]) - g(ν_ith_deg[end, "DATE"])
            Y2 = ν_ith_deg[end, "VALUE"] + rand(Normal(μ * Δft, σ * sqrt(Δgt)))
            Y1 = Y2 - ρ_u_ν_i * (Y2 - Y1)

            if ν[i] + 1 != ν[i+1]
                for l in ν[i]+1:ν[i+1]-1
                    u_l = maint[l, "TYPE"]#N'y-aurait-il pas un l+1 au lieu de l ?
                    ρ_u_l = type_to_ρ[u_l]
                    Δfτ_l = f(maint[l, "DATE"]) - f(maint[l - 1, "DATE"])
                    Δgτ_l = g(maint[l, "DATE"]) - g(maint[l - 1, "DATE"])
                    Y2 = Y1 + rand(Normal(μ * Δfτ_l, σ * sqrt(Δgτ_l)))
                    Y1 = Y2 - ρ_u_l * (Y2 - Y1)
                end
            end

            Δft = f(ν_iplus1th_first_deg[1, "DATE"]) - f(maint[ν[i], "DATE"])
            Δgt = g(ν_iplus1th_first_deg[1, "DATE"]) - g(maint[ν[i], "DATE"])
            ν_iplus1th_first_deg[1, "VALUE"] = Y1 + rand(Normal(μ * Δft, σ * sqrt(Δgt)))
            
            new_deg = vcat(new_deg, ν_ith_deg[2:end, "VALUE"], ν_iplus1th_first_deg.VALUE)
        else
            new_deg = vcat(new_deg, ν_ith_deg[2:end, "VALUE"])
        end
    end

    degradationdata.degradations.VALUE = new_deg

    return degradationdata
end

#####GENERATEUR DE DONNEES ALEATOIRES NON LINEAIRE INSPIRE DE LAURENT#####
function which_maint_type(type::String, maintTypes::Vector{String})
    if isempty(maintTypes)
        return 1
    else
        for (i,mtype) in enumerate(maintTypes)
            if mtype == type
                return i
            end
        end
    end
    @warn string("The maintenance type ", type, " apear in the data and has not been defined, consequantly it is assimilated to type ", maintTypes[1])
    return 1
end

function generateur_laurent_non_lin(model::WienerARD∞, data::DegradationData, from::Int64=0, f::Function=identity, g::Function=identity)
    try
        # Attempt to call `f` with Float64 input
        result = f(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Drift function `f` must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `f` failed, clarify the error
        throw(ArgumentError("Drift function `f` must accept Float64 input"))
    end

    try
        # Attempt to call `f` with Float64 input
        result = g(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Dispersion function `g` must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `f` failed, clarify the error
        throw(ArgumentError("Dispersion function `g` must accept Float64 input"))
    end

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
            dtf = [f(x) for x in t[2:end]] - [f(x) for x in t[1:(end-1)]]
            dtg = [g(x) for x in t[2:end]] - [g(x) for x in t[1:(end-1)]]
            rd = rand(Normal(), length(dtf))
            #res[i] = rd
            degval = deg0 .+ cumsum(rd .* modelParams[2] .* sqrt.(dtg) .+ modelParams[1] .* dtf)
            deg.VALUE[i_maint] = degval[1:(end-ip)]
            deg0 = (1 - ρ) * degval[end]
        end
    end

    data.infos["μ"] = modelParams[1]
    data.infos["σ"] = modelParams[2]
    data.infos["ρ"] = modelParams[3:end]
    data.infos["ρTypes"] = maintTypes
    data.infos["simulationModel"] = "ARDinf-Wiener"

    return data
end

function generateur_laurent_non_lin(model::WienerARD1, data::DegradationData, from::Int64=0, f::Function=identity, g::Function=identity)
    try
        # Attempt to call `f` with Float64 input
        result = f(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Drift function `f` must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `f` failed, clarify the error
        throw(ArgumentError("Drift function `f` must accept Float64 input"))
    end

    try
        # Attempt to call `f` with Float64 input
        result = g(1.)

        # Check if the output is Float64
        if !(result isa Float64)
            throw(ArgumentError("Dispersion function `g` must return Float64, got $(typeof(result))"))
        end
    catch e
        # If calling `f` failed, clarify the error
        throw(ArgumentError("Dispersion function `g` must accept Float64 input"))
    end

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
            dtf = [f(x) for x in t[2:end]] - [f(x) for x in t[1:(end-1)]]
            dtg = [g(x) for x in t[2:end]] - [g(x) for x in t[1:(end-1)]]
            rd = rand(Normal(), length(dtf))
            #res[i] = rd
            degval = deg0 .+ cumsum(rd .* modelParams[2] .* sqrt.(dtg) .+ modelParams[1] .* dtf)
            deg.VALUE[i_maint] = degval[1:(end-ip)]
            deg0 = degval[end] - ρ * (degval[end] - deg0)
        end
    end

    data.infos["μ"] = modelParams[1]
    data.infos["σ"] = modelParams[2]
    data.infos["ρ"] = modelParams[3:end]
    data.infos["ρTypes"] = maintTypes
    data.infos["simulationModel"] = "ARD1-Wiener"

    return data
end