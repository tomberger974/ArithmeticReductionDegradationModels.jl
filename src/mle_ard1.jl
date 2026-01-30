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