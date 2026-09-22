#####FUNCTIONS FOR SIMULATIONS#####
function tests(wienerARD1::WienerARD1, degradationdata::DegradationData, nb_tests::Int64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations; approprie pour les structures assymetriques
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)

    for i in 1:nb_tests
        rand!(wienerARD1, degradationdata)
        
        ith_wienerARD1 = fit_mle(wienerARD1, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion
end

function tests(wienerARD∞::WienerARD∞, degradationdata::DegradationData, nb_tests::Int64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD∞.maintenances), nb_tests)

    for i in 1:nb_tests
        rand!(wienerARD∞, degradationdata)
        
        ith_wienerARD∞ = fit_mle(wienerARD∞, degradationdata)

        mle_tests_maintenance_effect[:, i] = collect(ith_wienerARD∞[1, 3:end])
        mle_tests_drift[i] = ith_wienerARD∞[1, "μ"]
        mle_tests_dispersion[i] = ith_wienerARD∞[1, "σ"]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion
end

function tests(wienerARD1::WienerARD1, multiple_degradationdata::Vector{DegradationData}, nb_tests::Int64)#version où il suffit de rentrer une structure de degradationdata
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)

    for i in 1:nb_tests
        for degradationdata in multiple_degradationdata
            rand!(wienerARD1, degradationdata, 0)
        end

        ith_wienerARD1 = fit_mle(wienerARD1, multiple_degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion
end

function tests(wienerARD1::WienerARD1, k::Int64, nj::Int64, Δt::Float64, nb_tests::Int64)#version avec des maintenances et des observations periodiques
    T = (k+1)*nj
    τ = Vector{Float64}(undef, k)
    for i in 1:k
        τ[i] = i*(nj+1)
    end
    τ_types = ["P" for i in 1:k]
    T = τ[end] + nj
    degradationdata = CreateData2(τ, Δt, T, τ_types)
    filter!(row -> row.TYPE in ["Between"], degradationdata.degradations)
    
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)

    for i in 1:nb_tests
        rand!(wienerARD1, degradationdata, 0)
        
        ith_wienerARD1 = fit_mle(wienerARD1, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
    end
    
    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion
end



#####FUNCTIONS TO PLOT THE SIMULATIONS#####
function plot_tests(wienerARD∞::WienerARD∞, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion; savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD∞.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD∞.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    ylims!(-.5, 4.5)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    ylims!(0., 12.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)
    
    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\"*savename*".png", fig)
    end

    display(fig)
end