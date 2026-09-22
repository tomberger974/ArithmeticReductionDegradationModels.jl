#####MONO-SYSTEM POWER FUNCTION ESTIMATORS ARD∞#####
function coefficient(ρ::Vector{Float64}, degradationdata::DegradationData, used_wienerARD∞::WienerARD∞, α::Float64, β::Float64)
    deg = vcat(DataFrame(DATE = 0., NB_MAINTENANCES = 0, VALUE = 0., TYPE = "Between"), degradationdata.degradations)#We create a fictive observation at time 0 with no effect to simplify the code
    maint = degradationdata.maintenances
    
    type_to_ρ = Dict(used_wienerARD∞.maintenances.TYPE .=> ρ)
    
    NBM = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint))
    ν = NBM .+ 1
    K = length(ν)
    
    Y = Vector{Float64}(undef, K - 1)
    E = Vector{Float64}(undef, K - 1)
    V = Vector{Float64}(undef, K - 1)
    
    for i in 1:K-1
        u_ν = [maint[l, "TYPE"] for l in ν[i]:ν[i+1]-1]
        P = prod(1 .- [type_to_ρ[t] for t in u_ν])
        
        
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i] - 1, deg)#keep only the observations between \nu[i]-1 and \nu[i]
        ν_iplus1th_first_deg = deg[findfirst(deg.NB_MAINTENANCES .== ν[i + 1] - 1), :]#find the first observation following all the observations stored in ν_ith_deg
        
        Y[i] = ν_iplus1th_first_deg.VALUE - ν_ith_deg.VALUE[end] * P
        
        Δt_ν_iplus1_1_α = ν_iplus1th_first_deg.DATE ^ α - maint[ν[i + 1] - 1, "DATE"] ^ α
        Δt_ν_iplus1_1_β = ν_iplus1th_first_deg.DATE ^ β - maint[ν[i + 1] - 1, "DATE"] ^ β
        Δt_ν_i_nνiplus1_α = maint[ν[i], "DATE"] ^ α - ν_ith_deg[end, "DATE"] ^ α
        Δt_ν_i_nνiplus1_β = maint[ν[i], "DATE"] ^ β - ν_ith_deg[end, "DATE"] ^ β
        E[i] = Δt_ν_i_nνiplus1_α * P + Δt_ν_iplus1_1_α
        V[i] = Δt_ν_i_nνiplus1_β * P^2 + Δt_ν_iplus1_1_β
        
        if ν[i+1] != ν[i] + 1
            for l in ν[i]+1:ν[i+1]-1
                u_l = [maint[m, "TYPE"] for m in l:ν[i+1]-1]
                P2 = prod(1 .- [type_to_ρ[t] for t in u_l])
                Δτ_l_α = maint[l, "DATE"] ^ α - maint[l - 1, "DATE"] ^ α
                Δτ_l_β = maint[l, "DATE"] ^ β - maint[l - 1, "DATE"] ^ β
                E[i] += P2 * Δτ_l_α
                V[i] += P2^2 * Δτ_l_β
            end
        end
    end
    
    return [Y, E, V]
end

function fit_mle_drift∞(degradationdata::DegradationData, Y::Vector{Float64}, E::Vector{Float64}, V::Vector{Float64}, α::Float64, β::Float64)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    
    ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
    K = length(ν)
    
    S1, S2 = 0., 0.
    
    if ν[1] == 0#if there is an observation before first maintenance, need to add first_obs - 0 as increment
        S1 += deg[1, "VALUE"] * deg[1, "DATE"] ^ (α - β)
        S2 += deg[1, "DATE"] ^ (2 * α - β)
    end
    
    for i in 1:K
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
        if nrow(ν_ith_deg) != 1
            ΔY_ν_i = diff(ν_ith_deg.VALUE)
            Δt_ν_i_α = diff(ν_ith_deg.DATE .^ α)
            Δt_ν_i_β = diff(ν_ith_deg.DATE .^ β)
            S1 += sum(ΔY_ν_i .* Δt_ν_i_α ./ Δt_ν_i_β)
            S2 += sum(Δt_ν_i_α .^ 2 ./ Δt_ν_i_β)
        end
    end
    
    return (S1 + sum(E .* Y ./ V)) / (S2 + sum(E .^ 2 ./ V))
end

function fit_mle_dispersion∞(degradationdata::DegradationData, Y::Vector{Float64}, E::Vector{Float64}, V::Vector{Float64}, drift::Float64, α::Float64, β::Float64)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    
    ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
    K = length(ν)
    N = nrow(deg)
    
    S = 0.
    if ν[1] == 0
        S += (deg[1, "VALUE"] - drift * deg[1, "DATE"] ^ α)^2 / deg[1, "DATE"] ^ β
    end
    
    for i in 1:K
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
        if nrow(ν_ith_deg) != 1
            ΔY_ν_i = diff(ν_ith_deg.VALUE)
            Δt_ν_i_α = diff(ν_ith_deg.DATE .^ α)
            Δt_ν_i_β = diff(ν_ith_deg.DATE .^ β)
            S += sum((ΔY_ν_i .- drift .* Δt_ν_i_α).^2 ./ Δt_ν_i_β)
        end
    end
    
    return (S + sum((Y .- drift .* E) .^ 2 ./ V)) / N
end

function fit_mle_maintenance_effect_α_β(used_wienerARD∞::WienerARD∞, degradationdata::DegradationData)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    N = nrow(deg)

    # Define the bounds
    U = length(used_wienerARD∞.maintenances.VALUE)
    lower = vcat([-Inf for i in 1:U], [0.0, 0.0])
    upper = [Inf for i in 1:U + 2]

    initial_ρ_α_β = vcat(used_wienerARD∞.maintenances.VALUE, 1, 1)
    
    function f(ρ_α_β)
        ρ = ρ_α_β[1:end-2]
        α = ρ_α_β[end-1]
        β = ρ_α_β[end]
        Y, E, V = coefficient(ρ, degradationdata, used_wienerARD∞, α, β)[1:3]
        
        drift = fit_mle_drift∞(degradationdata, Y, E, V, α, β)
        dispersion = fit_mle_dispersion∞(degradationdata, Y, E, V, drift, α, β)
        
        ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
        K = length(ν)
        
        S = 0.
        if ν[1] == 0
            S += log(deg[1, "DATE"] ^ β)
        end
        
        for i in 1:K
            ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
            if nrow(ν_ith_deg) != 1
                Δt_ν_i_β = diff(ν_ith_deg.DATE .^ β)
                S += sum((log.(Δt_ν_i_β)))
            end
        end
        
        return N * log(dispersion) + sum(log.(V)) + S
    end
    
    return optimize(f, lower, upper, initial_ρ_α_β, Fminbox())
end

function fit_mle_α_β(wienerARD∞::WienerARD∞, degradationdata::DegradationData)
    used_maintenance_indices = [type in degradationdata.maintenances.TYPE for type in wienerARD∞.maintenances.TYPE]
    used_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ_α_β = fit_mle_maintenance_effect_α_β(used_wienerARD∞, degradationdata).minimizer
    ρ = ρ_α_β[1:end-2]
    α = ρ_α_β[end-1]
    β = ρ_α_β[end]
    
    Y, E, V = coefficient(ρ, degradationdata, used_wienerARD∞, α, β)[1:3]
    
    drift = fit_mle_drift∞(degradationdata, Y, E, V, α, β)
    dispersion = sqrt(fit_mle_dispersion∞(degradationdata, Y, E, V, drift, α, β))
    
    new_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances.TYPE, drift, dispersion)
    new_wienerARD∞.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD∞.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD∞, α, β
end

function fit_mle_maintenance_effect_α(used_wienerARD∞::WienerARD∞, degradationdata::DegradationData)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    N = nrow(deg)
    initial_ρ_α = vcat(used_wienerARD∞.maintenances.VALUE, 1)

    U = length(used_wienerARD∞.maintenances.VALUE)
    lower = vcat([-Inf for i in 1:U], [0.0])
    upper = [Inf for i in 1:U + 1]
    
    function f(ρ_α)
        ρ = ρ_α[1:end-1]
        α = ρ_α[end]
        Y, E, V = coefficient(ρ, degradationdata, used_wienerARD∞, α, 1.)[1:3]
        
        drift = fit_mle_drift∞(degradationdata, Y, E, V, α, 1.)
        dispersion = fit_mle_dispersion∞(degradationdata, Y, E, V, drift, α, 1.)
        
        return N * log(dispersion) + sum(log.(V))
    end
    
    return optimize(f, lower, upper, initial_ρ_α, Fminbox())
end

function fit_mle_α(wienerARD∞::WienerARD∞, degradationdata::DegradationData)
    used_maintenance_indices = [type in degradationdata.maintenances.TYPE for type in wienerARD∞.maintenances.TYPE]
    used_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ_α = fit_mle_maintenance_effect_α(used_wienerARD∞, degradationdata).minimizer
    ρ = ρ_α[1:end-1]
    α = ρ_α[end]
    
    Y, V, E = coefficient(ρ, degradationdata, used_wienerARD∞, α, 1.)[1:3]
    
    drift = fit_mle_drift∞(degradationdata, Y, V, E, α, 1.)
    dispersion = sqrt(fit_mle_dispersion∞(degradationdata, Y, V, E, drift, α, 1.))
    
    new_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances.TYPE, drift, dispersion)
    new_wienerARD∞.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD∞.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD∞, α
end

function fit_mle_maintenance_effect_β(used_wienerARD∞::WienerARD∞, degradationdata::DegradationData)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    N = nrow(deg)

    # Define the bounds
    U = length(used_wienerARD∞.maintenances.VALUE)
    lower = vcat([-Inf for i in 1:U], [0.0])
    upper = [Inf for i in 1:U+1]

    initial_ρ_β = vcat(used_wienerARD∞.maintenances.VALUE, 1)
    
    function f(ρ_β)
        ρ = ρ_β[1:end-1]
        β = ρ_α_β[end]
        Y, E, V = coefficient(ρ, degradationdata, used_wienerARD∞, 1., β)[1:3]
        
        drift = fit_mle_drift∞(degradationdata, Y, E, V, 1., β)
        dispersion = fit_mle_dispersion∞(degradationdata, Y, E, V, drift, 1., β)
        
        ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
        K = length(ν)
        
        S = 0.
        if ν[1] == 0
            S += log(deg[1, "DATE"] ^ β)
        end
        
        for i in 1:K
            ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
            if nrow(ν_ith_deg) != 1
                Δt_ν_i_β = diff(ν_ith_deg.DATE .^ β)
                S += sum((log.(Δt_ν_i_β)))
            end
        end
        
        return N * log(dispersion) + sum(log.(V)) + S
    end
    
    return optimize(f, lower, upper, initial_ρ_β, Fminbox())
end

function fit_mle_β(wienerARD∞::WienerARD∞, degradationdata::DegradationData)
    used_maintenance_indices = [type in degradationdata.maintenances.TYPE for type in wienerARD∞.maintenances.TYPE]
    used_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ_β = fit_mle_maintenance_effect_β(used_wienerARD∞, degradationdata).minimizer
    ρ = ρ_β[1:end-1]
    β = ρ_β[end]
    
    Y, E, V = coefficient(ρ, degradationdata, used_wienerARD∞, 1., β)[1:3]
    
    drift = fit_mle_drift∞(degradationdata, Y, E, V, 1., β)
    dispersion = sqrt(fit_mle_dispersion∞(degradationdata, Y, E, V, drift, 1., β))
    
    new_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances.TYPE, drift, dispersion)
    new_wienerARD∞.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD∞.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD∞, β
end



#####MULTI-SYSTEM POWER FUNCTION ESTIMATORS ARD∞#####
function coefficient(ρ::Vector{Float64}, multiple_degradationdata::Vector{DegradationData}, used_wienerARD∞::WienerARD∞, α::Float64, β::Float64)
    nb_systeme = length(multiple_degradationdata)
    multiple_Y = Vector{Vector{Float64}}(undef, nb_systeme)
    multiple_E = Vector{Vector{Float64}}(undef, nb_systeme)
    multiple_V = Vector{Vector{Float64}}(undef, nb_systeme)
    
    for i in eachindex(multiple_degradationdata)
        multiple_Y[i], multiple_E[i], multiple_V[i] = coefficient(ρ, multiple_degradationdata[i], used_wienerARD∞, α, β)[1:3]
    end
    
    return [multiple_Y, multiple_E, multiple_V]
end

function fit_mle_drift∞(multiple_degradationdata::Vector{DegradationData}, multiple_Y::Vector{Vector{Float64}}, multiple_E::Vector{Vector{Float64}}, multiple_V::Vector{Vector{Float64}}, α::Float64, β::Float64)
    S1, S2 = 0., 0.
    
    for l in eachindex(multiple_degradationdata)
        deg = multiple_degradationdata[l].degradations
        maint = multiple_degradationdata[l].maintenances
        
        ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
        K = length(ν)
        
        if ν[1] == 0#if there is an observation before first maintenance, need to add first_obs - 0 as increment
            S1 += deg.VALUE[1] * deg.DATE[1] ^ (α - β)
            S2 += deg.DATE[1] ^ (2 * α - β)
        end
        
        for i in 1:K
            ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
            if nrow(ν_ith_deg) != 1
                ΔY_ν_i = diff(ν_ith_deg.VALUE)
                Δt_ν_i_α = diff(ν_ith_deg.DATE .^ α)
                Δt_ν_i_β = diff(ν_ith_deg.DATE .^ β)
                S1 += sum(ΔY_ν_i .* Δt_ν_i_α ./ Δt_ν_i_β)
                S2 += sum(Δt_ν_i_α .^ 2 ./ Δt_ν_i_β)
            end
        end

        S1 += sum(multiple_E[l] .* multiple_Y[l] ./ multiple_V[l])
        S2 + sum(multiple_E[l] .^2 ./multiple_V[l])
    end
    
    return S1 / S2
end

function fit_mle_dispersion∞(multiple_degradationdata::Vector{DegradationData}, multiple_Y::Vector{Vector{Float64}}, multiple_E::Vector{Vector{Float64}}, multiple_V::Vector{Vector{Float64}}, drift::Float64, α::Float64, β::Float64)
    S = 0.
    N = 0
    
    for l in eachindex(multiple_degradationdata)
        deg = multiple_degradationdata[l].degradations
        maint = multiple_degradationdata[l].maintenances
        N += nrow(deg)
        
        ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
        K = length(ν)
        
        if ν[1] == 0
            S += (deg.VALUE[1] - drift * deg.DATE[1]^α)^2 / deg.DATE[1]^β
        end
        
        for i in 1:K
            ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
            if nrow(ν_ith_deg) != 1
                Δt_ν_i_α = diff(ν_ith_deg.DATE .^ α)
                Δt_ν_i_β = diff(ν_ith_deg.DATE .^ β)
                ΔY_ν_i = diff(ν_ith_deg.VALUE)
                S += sum((ΔY_ν_i .- drift .* Δt_ν_i_α) .^ 2 ./ Δt_ν_i_β)
            end
        end
        
        S += sum((multiple_Y[l] .- drift .* multiple_E[l]) .^ 2 ./ multiple_V[l])
    end

    return S / N
end

function fit_mle_maintenance_effect_α_β(used_wienerARD∞::WienerARD∞, multiple_degradationdata::Vector{DegradationData})
    N = 0
    for degradationdata in multiple_degradationdata
        N += nrow(degradationdata.degradations)
    end
    
    initial_ρ_α_β = vcat(used_wienerARD∞.maintenances.VALUE, 1., 1.)

    U = length(used_wienerARD∞.maintenances.VALUE)
    lower = vcat([-Inf for i in 1:U], [0., 0.])
    upper = [Inf for i in 1:U + 2]
    
    function f(ρ_α_β)
        ρ = ρ_α_β[1:end-2]
        α = ρ_α_β[end-1]
        β = ρ_α_β[end]
        multiple_Y, multiple_E, multiple_V = coefficient(ρ, multiple_degradationdata, used_wienerARD∞, α, β)[1:3]
        
        drift = fit_mle_drift∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, α, β)
        dispersion = fit_mle_dispersion∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, drift, α, β)
        
        return N * log(dispersion) + sum([sum(log.(V)) for V in multiple_V]) + sum([sum_Δ_t_β(degradationdata, β) for degradationdata in multiple_degradationdata])
    end
    
    return optimize(f, lower, upper, initial_ρ_α_β, Fminbox())
end

function fit_mle_α_β(wienerARD∞::WienerARD∞, multiple_degradationdata::Vector{DegradationData})
    used_maintenance_indices = [type in reduce(vcat, [degradationdata.maintenances.TYPE for degradationdata in multiple_degradationdata]) for type in wienerARD∞.maintenances.TYPE]
    used_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ_α_β = fit_mle_maintenance_effect_α_β(used_wienerARD∞, multiple_degradationdata).minimizer
    ρ = ρ_α_β[1:end-2]
    α = ρ_α_β[end-1]
    β = ρ_α_β[end]
    
    multiple_Y, multiple_E, multiple_V = coefficient(ρ, multiple_degradationdata, used_wienerARD∞, α, β)[1:3]
    drift = fit_mle_drift∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, α, β)
    dispersion = sqrt(fit_mle_dispersion∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, drift, α, β))
    
    new_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances.TYPE, drift, dispersion)
    new_wienerARD∞.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD∞.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD∞, α, β
end

function fit_mle_maintenance_effect_α(used_wienerARD∞::WienerARD∞, multiple_degradationdata::Vector{DegradationData})
    N = 0
    for degradationdata in multiple_degradationdata
        N += nrow(degradationdata.degradations)
    end
    
    initial_ρ_α = vcat(used_wienerARD∞.maintenances.VALUE, 1.)

    U = length(used_wienerARD∞.maintenances.VALUE)
    lower = vcat([-Inf for i in 1:U], [0.0])
    upper = [Inf for i in 1:U + 1]
    
    function f(ρ_α)
        ρ = ρ_α[1:end-1]
        α = ρ_α[end]
        multiple_Y, multiple_E, multiple_V = coefficient(ρ, multiple_degradationdata, used_wienerARD∞, α, 1.)[1:3]
        
        drift = fit_mle_drift∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, α, 1.)
        dispersion = fit_mle_dispersion∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, drift, α, 1.)

        return N * log(dispersion) + sum([sum(log.(V)) for V in multiple_V])
    end
    
    return optimize(f, lower, upper, initial_ρ_α, Fminbox())
end

function fit_mle_α(wienerARD∞::WienerARD∞, multiple_degradationdata::Vector{DegradationData})
    used_maintenance_indices = [type in reduce(vcat, [degradationdata.maintenances.TYPE for degradationdata in multiple_degradationdata]) for type in wienerARD∞.maintenances.TYPE]
    used_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ_α = fit_mle_maintenance_effect_α(used_wienerARD∞, multiple_degradationdata).minimizer
    ρ = ρ_α[1:end-1]
    α = ρ_α[end]
    
    multiple_Y, multiple_E, multiple_V = coefficient(ρ, multiple_degradationdata, used_wienerARD∞, α, 1.)[1:3]
    drift = fit_mle_drift∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, α, 1.)
    dispersion = sqrt(fit_mle_dispersion∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, drift, α, 1.))
    
    new_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances.TYPE, drift, dispersion)
    new_wienerARD∞.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD∞.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD∞, α
end

function fit_mle_maintenance_effect_β(used_wienerARD∞::WienerARD∞, multiple_degradationdata::Vector{DegradationData})
    N = 0
    for degradationdata in multiple_degradationdata
        N += nrow(degradationdata.degradations)
    end
    
    initial_ρ_β = vcat(used_wienerARD∞.maintenances.VALUE, 1.)

    U = length(used_wienerARD∞.maintenances.VALUE)
    lower = vcat([-Inf for i in 1:U], [0.])
    upper = [Inf for i in 1:U+1]
    
    function f(ρ_β)
        ρ = ρ_β[1:end-1]
        β = ρ_β[end]
        multiple_Y, multiple_E, multiple_V = coefficient(ρ, multiple_degradationdata, used_wienerARD∞, 1., β)[1:3]

        drift = fit_mle_drift∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, 1., β)
        dispersion = fit_mle_dispersion∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, drift, 1., β)
        
        return N * log(dispersion) + sum([sum(log.(V)) for V in multiple_V]) + sum([sum_Δ_t_β(degradationdata, β) for degradationdata in multiple_degradationdata])
    end
    
    return optimize(f, lower, upper, initial_ρ_β, Fminbox())
end

function fit_mle_β(wienerARD∞::WienerARD∞, multiple_degradationdata::Vector{DegradationData})
    used_maintenance_indices = [type in reduce(vcat, [degradationdata.maintenances.TYPE for degradationdata in multiple_degradationdata]) for type in wienerARD∞.maintenances.TYPE]
    used_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ_β = fit_mle_maintenance_effect_β(used_wienerARD∞, multiple_degradationdata).minimizer
    ρ = ρ_β[1:end-1]
    β = ρ_β[end]
    
    multiple_Y, multiple_E, multiple_V = coefficient(ρ, multiple_degradationdata, used_wienerARD∞, 1., β)[1:3]
    drift = fit_mle_drift∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, 1., β)
    dispersion = sqrt(fit_mle_dispersion∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, drift, 1., β))
    
    new_wienerARD∞ = WienerARD∞(wienerARD∞.maintenances.TYPE, drift, dispersion)
    new_wienerARD∞.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD∞.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD∞, β
end

function sum_Δ_t_β(degradationdata::DegradationData, β::Float64)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances
    ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
    K = length(ν)
    
    S = 0.
    
    if ν[1] == 0
        S += log(deg[1, "DATE"] ^ β)
    end
    
    for i in 1:K
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
        if nrow(ν_ith_deg) != 1
            Δt_ν_i_β = diff(ν_ith_deg.DATE .^ β)
            S += sum((log.(Δt_ν_i_β)))
        end
    end

    return S
end