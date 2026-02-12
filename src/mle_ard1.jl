#####MONO-SYSTEM ESTIMATORS#####
function coefficient(ρ::Vector{Float64}, degradationdata::DegradationData, used_wienerARD1::WienerARD1)
    deg = vcat(DataFrame(DATE = 0., NB_MAINTENANCES = 0, VALUE = 0., TYPE = "Between"), degradationdata.degradations)
    maint = degradationdata.maintenances
    
    type_to_ρ = Dict(used_wienerARD1.maintenances.TYPE .=> ρ)
    
    NBM = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint))
    ν = NBM .+ 1
    K = length(ν)
    
    a = Vector{Float64}(undef, K - 1)
    M = Vector{Float64}(undef, K - 1)
    V = Vector{Float64}(undef, K - 1)
    C = Vector{Float64}(undef, K - 2)
    
    for i in 1:K-1
        u_ν_i = maint.TYPE[ν[i]]
        ρ_u_ν_i = type_to_ρ[u_ν_i]
        
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i] - 1, deg)#keep only the observations between \nu[i]-1 and \nu[i]
        ν_iplus1th_first_deg = deg[findfirst(deg.NB_MAINTENANCES .== ν[i + 1] - 1), :]#find the first observation following all the observations stored in ν_ith_deg
        
        z_ν_i = ν_iplus1th_first_deg.VALUE - ν_ith_deg.VALUE[end]#computes the observed jump
        a[i] = nrow(ν_ith_deg) > 1 ? z_ν_i + ρ_u_ν_i * sum(diff(ν_ith_deg.VALUE)) : z_ν_i#computes the a associated to a given z
        
        Δt_ν_iplus1 = ν_iplus1th_first_deg.DATE - maint.DATE[ν[i + 1] - 1]
        Δt_ν_i_end = maint.DATE[ν[i]] - ν_ith_deg.DATE[end]
        M[i] = (1 - ρ_u_ν_i) * Δt_ν_i_end + Δt_ν_iplus1
        V[i] = (1 - ρ_u_ν_i)^2 * Δt_ν_i_end + Δt_ν_iplus1
        if ν[i] != 1
            Δt_ν_i_1 = ν_ith_deg.DATE[1] - maint.DATE[ν[i] - 1]
            M[i] += - ρ_u_ν_i * Δt_ν_i_1
            V[i] += ρ_u_ν_i^2 * Δt_ν_i_1
        end
        
        if i != K - 1#there is one less C defined
            u_ν_iplus1 = maint.TYPE[ν[i+1]]
            ρ_u_ν_iplus1 = type_to_ρ[u_ν_iplus1]
            C[i] = ρ_u_ν_iplus1 *  Δt_ν_iplus1
        end
        
        if ν[i+1] != ν[i] + 1
            for l in ν[i]+1:ν[i+1]-1
                u_l = maint.TYPE[l]#N'y-aurait-il pas un l+1 au lieu de l ?
                ρ_u_l = type_to_ρ[u_l]
                Δτ_l = maint.DATE[l] - maint.DATE[l - 1]
                M[i] += (1 - ρ_u_l) * Δτ_l
                V[i] += (1 - ρ_u_l)^2 * Δτ_l
            end
        end
    end
    
    return [a, M, V, C]
end

function fit_mle_drift(degradationdata::DegradationData, a::Vector{Float64}, M::Vector{Float64}, Σ_minus1::Matrix{Float64})
    deg = vcat(DataFrame(DATE = 0., NB_MAINTENANCES = 0, VALUE = 0., TYPE = "Between"), degradationdata.degradations)
    maint = degradationdata.maintenances
    
    ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
    K = length(ν)
    
    S1, S2 = 0., 0.
    
    for i in 1:K
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
        if nrow(ν_ith_deg) != 1
            S1 += ν_ith_deg.VALUE[end] - ν_ith_deg.VALUE[1]
            S2 += ν_ith_deg.DATE[end] - ν_ith_deg.DATE[1]
        end
    end
    
    return (S1 + dot(M, Σ_minus1 * a)) / (S2 + dot(M, Σ_minus1 * M))
end

function fit_mle_dispersion(degradationdata::DegradationData, a::Vector{Float64}, M::Vector{Float64}, Σ_minus1::Matrix{Float64}, drift::Float64)
    deg = vcat(DataFrame(DATE = 0., NB_MAINTENANCES = 0, VALUE = 0., TYPE = "Between"), degradationdata.degradations)
    maint = degradationdata.maintenances
    
    NBM = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint))
    ν = NBM .+ 1
    K = length(ν)
    N = nrow(deg) - 1#we added a virtual observation at date 0. so we need to remove it in the counting
    
    S = 0.
    
    for i in 1:K
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i] - 1, deg)#keep only the observations between \nu[i]-1 and \nu[i]
        if nrow(ν_ith_deg) != 1
            Δt_ν_i = diff(ν_ith_deg.DATE)
            ΔY_ν_i = diff(ν_ith_deg.VALUE)
            S += sum((ΔY_ν_i .- (drift .* Δt_ν_i)).^2 ./ Δt_ν_i)
        end
    end
    
    vecteur = a .- drift .* M
    S += dot(vecteur, Σ_minus1 * vecteur)
    
    return S / N
end

function fit_mle_maintenance_effect(used_wienerARD1::WienerARD1, degradationdata::DegradationData)
    N = nrow(degradationdata.degradations)
    initial_ρ = used_wienerARD1.maintenances.VALUE
    
    function f(ρ)
        a, M, V, C = coefficient(ρ, degradationdata, used_wienerARD1)[1:4]
        Σ = Tridiagonal(-C, V, -C)
        Σ_minus1 = inv(Σ)
        
        drift = fit_mle_drift(degradationdata, a, M, Σ_minus1)
        dispersion = fit_mle_dispersion(degradationdata, a, M, Σ_minus1, drift)
        
        return N * log(dispersion) + log(det(Σ))
    end
    
    return optimize(f, initial_ρ)
end

function fit_mle(wienerARD1::WienerARD1, degradationdata::DegradationData)
    used_maintenance_indices = [type in degradationdata.maintenances.TYPE for type in wienerARD1.maintenances.TYPE]
    used_wienerARD1 = WienerARD1(wienerARD1.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ = fit_mle_maintenance_effect(used_wienerARD1, degradationdata).minimizer
    
    a, M, V, C = coefficient(ρ, degradationdata, used_wienerARD1)[1:4]
    Σ = Tridiagonal(-C, V, -C)
    Σ_minus1 = inv(Σ)
    
    drift = fit_mle_drift(degradationdata, a, M, Σ_minus1)
    dispersion = sqrt(fit_mle_dispersion(degradationdata, a, M, Σ_minus1, drift))
    
    new_wienerARD1 = WienerARD1(wienerARD1.maintenances.TYPE, drift, dispersion)
    new_wienerARD1.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD1.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD1
end

function fit_mle!(wienerARD1::WienerARD1, degradationdata::DegradationData)
    wienerARD1 = fit_mle(wienerARD1, degradationdata)
    
    return wienerARD1
end

#####MULTI-SYSTEM ESTIMATORS#####
function coefficient(ρ::Vector{Float64}, multiple_degradationdata::Vector{DegradationData}, used_wienerARD1::WienerARD1)
    nb_systeme = length(multiple_degradationdata)
    multiple_a = Vector{Vector{Float64}}(undef, nb_systeme)
    multiple_M = Vector{Vector{Float64}}(undef, nb_systeme)
    multiple_V = Vector{Vector{Float64}}(undef, nb_systeme)
    multiple_C = Vector{Vector{Float64}}(undef, nb_systeme)
    
    for i in eachindex(multiple_degradationdata)
        multiple_a[i], multiple_M[i], multiple_V[i], multiple_C[i] = coefficient(ρ, multiple_degradationdata[i], used_wienerARD1)[1:4]
    end
    
    return [multiple_a, multiple_M, multiple_V, multiple_C]
end

function fit_mle_drift(multiple_degradationdata::Vector{DegradationData}, multiple_a::Vector{Vector{Float64}}, multiple_M::Vector{Vector{Float64}}, multiple_Σ_minus1::Vector{Matrix{Float64}})
    S1, S2 = 0., 0.
    
    for l in eachindex(multiple_degradationdata)
        deg = multiple_degradationdata[l].degradations
        maint = multiple_degradationdata[l].maintenances
        
        ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
        K = length(ν)
        
        if ν[1] == 0#if there is an observation before first maintenance, need to add first_obs - 0 as increment
            S1 += deg.VALUE[1]
            S2 += deg.DATE[1]
        end
        
        for i in 1:K
            ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
            if nrow(ν_ith_deg) != 1
                S1 += ν_ith_deg.VALUE[end] - ν_ith_deg.VALUE[1]
                S2 += ν_ith_deg.DATE[end] - ν_ith_deg.DATE[1]
            end
        end
    end
    
    return (S1 + sum(dot.(multiple_M, multiple_Σ_minus1 .* multiple_a))) / (S2 + sum(dot.(multiple_M, multiple_Σ_minus1 .* multiple_M)))
end

function fit_mle_dispersion(multiple_degradationdata::Vector{DegradationData}, multiple_a::Vector{Vector{Float64}}, multiple_M::Vector{Vector{Float64}}, multiple_Σ_minus1::Vector{Matrix{Float64}}, drift::Float64)
    S = 0.
    N = 0
    
    for l in eachindex(multiple_degradationdata)
        deg = multiple_degradationdata[l].degradations
        maint = multiple_degradationdata[l].maintenances
        N += nrow(deg)
        
        ν = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
        K = length(ν)
        
        if ν[1] == 0
            S += (deg.VALUE[1] - drift * deg.DATE[1])^2 / deg.DATE[1]
        end
        
        for i in 1:K
            ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i], deg)#keep only the observations between \nu[i]-1 and \nu[i]
            if nrow(ν_ith_deg) != 1
                Δt_ν_i = diff(ν_ith_deg.DATE)
                ΔY_ν_i = diff(ν_ith_deg.VALUE)
                S += sum((ΔY_ν_i .- drift .* Δt_ν_i).^2 ./ Δt_ν_i)
            end
        end
        
        V = multiple_a[l] .- drift .* multiple_M[l]
        S += dot(V, multiple_Σ_minus1[l] * V)
    end
    
    return S / N
end

function fit_mle_maintenance_effect(used_wienerARD1::WienerARD1, multiple_degradationdata::Vector{DegradationData})
    N = 0
    for degradationdata in multiple_degradationdata
        N += nrow(degradationdata.degradations)
    end
    
    initial_ρ = used_wienerARD1.maintenances.VALUE
    
    function f(ρ)
        multiple_a, multiple_M, multiple_V, multiple_C = coefficient(ρ, multiple_degradationdata, used_wienerARD1)[1:4]
        multiple_Σ = Tridiagonal.(.-multiple_C, multiple_V, .-multiple_C)
        multiple_Σ_minus1 = inv.(multiple_Σ)
        
        drift = fit_mle_drift(multiple_degradationdata, multiple_a, multiple_M, multiple_Σ_minus1)
        dispersion = fit_mle_dispersion(multiple_degradationdata, multiple_a, multiple_M, multiple_Σ_minus1, drift)
        
        return N * log(dispersion) + sum(log.(det.(multiple_Σ)))
    end
    
    return optimize(f, initial_ρ)
end

function fit_mle(wienerARD1::WienerARD1, multiple_degradationdata::Vector{DegradationData})
    used_maintenance_indices = [type in reduce(vcat, [degradationdata.maintenances.TYPE for degradationdata in multiple_degradationdata]) for type in wienerARD1.maintenances.TYPE]
    used_wienerARD1 = WienerARD1(wienerARD1.maintenances[used_maintenance_indices, "TYPE"], 0., 0.)
    
    ρ = fit_mle_maintenance_effect(used_wienerARD1, multiple_degradationdata).minimizer
    
    multiple_a, multiple_M, multiple_V, multiple_C = coefficient(ρ, multiple_degradationdata, used_wienerARD1)[1:4]
    multiple_Σ = Tridiagonal.(.-multiple_C, multiple_V, .-multiple_C)
    multiple_Σ_minus1 = inv.(multiple_Σ)
    
    drift = fit_mle_drift(multiple_degradationdata, multiple_a, multiple_M, multiple_Σ_minus1)
    dispersion = sqrt(fit_mle_dispersion(multiple_degradationdata, multiple_a, multiple_M, multiple_Σ_minus1, drift))
    
    new_wienerARD1 = WienerARD1(wienerARD1.maintenances.TYPE, drift, dispersion)
    new_wienerARD1.maintenances.VALUE[used_maintenance_indices] = ρ
    new_wienerARD1.maintenances.VALUE[Not(used_maintenance_indices)] .= 0.
    
    return new_wienerARD1
end

function fit_mle!(wienerARD1::WienerARD1, multiple_degradationdata::Vector{DegradationData})
    wienerARD1 = fit_mle(wienerARD1, multiple_degradationdata)
    
    return wienerARD1
end