#####FUNCTIONS FOR SIMULATIONS#####
function tests_α_β(wienerARD1::WienerARD1, degradationdata::DegradationData, nb_tests::Int64, α::Float64, β::Float64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations; approprie pour les structures assymetriques
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)
    mle_tests_α = Vector{Float64}(undef, nb_tests)
    mle_tests_β = Vector{Float64}(undef, nb_tests)

    for i in 1:nb_tests
        generateur_laurent_non_lin(wienerARD1, degradationdata, 0, x -> x ^ α, x -> x ^ β)
        
        ith_wienerARD1, mle_tests_α[i], mle_tests_β[i] = fit_mle_α_β(wienerARD1, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α, mle_tests_β
end

function tests_α(wienerARD1::WienerARD1, degradationdata::DegradationData, nb_tests::Int64, α::Float64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations; approprie pour les structures assymetriques
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)
    mle_tests_α = Vector{Float64}(undef, nb_tests)

    for i in 1:nb_tests
        generateur_laurent_non_lin(wienerARD1, degradationdata, 0, x -> x ^ α, x -> x)
        
        ith_wienerARD1, mle_tests_α[i]= ARD.fit_mle_α(wienerARD1, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α
end

function tests_α(wienerARD1::WienerARD1, multiple_degradationdata::Vector{DegradationData}, nb_tests::Int64, α::Float64)
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)
    mle_tests_α = Vector{Float64}(undef, nb_tests)

    for i in 1:nb_tests
        for degradationdata in multiple_degradationdata
            generateur_laurent_non_lin(wienerARD1, degradationdata, 0, x -> x ^ α, x -> x)
        end
        
        ith_wienerARD1, mle_tests_α[i]= ARD.fit_mle_α(wienerARD1, multiple_degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α
end

function tests_β(wienerARD1::WienerARD1, degradationdata::DegradationData, nb_tests::Int64, β::Float64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations; approprie pour les structures assymetriques
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)
    mle_tests_β = Vector{Float64}(undef, nb_tests)

    for i in 1:nb_tests
        generateur_laurent_non_lin(wienerARD1, degradationdata, 0, x -> x, x -> x ^ β)
        
        ith_wienerARD1, mle_tests_β[i] = fit_mle_β(wienerARD1, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_β
end

function tests_α_β(wienerARD∞::WienerARD∞, degradationdata::DegradationData, nb_tests::Int64, α::Float64, β::Float64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations; approprie pour les structures assymetriques
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD∞.maintenances), nb_tests)
    mle_tests_α = Vector{Float64}(undef, nb_tests)
    mle_tests_β = Vector{Float64}(undef, nb_tests)

    for i in 1:nb_tests
        generateur_laurent_non_lin(wienerARD∞, degradationdata, 0, x -> x ^ α, x -> x ^ β)
        
        ith_wienerARD∞, mle_tests_α[i], mle_tests_β[i] = fit_mle_α_β(wienerARD∞, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD∞.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD∞.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD∞.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α, mle_tests_β
end

function tests_α(wienerARD∞::WienerARD∞, degradationdata::DegradationData, nb_tests::Int64, α::Float64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations; approprie pour les structures assymetriques
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD∞.maintenances), nb_tests)
    mle_tests_α = Vector{Float64}(undef, nb_tests)

    for i in 1:nb_tests
        generateur_laurent_non_lin(wienerARD∞, degradationdata, 0, x -> x ^ α, x -> x)
        
        ith_wienerARD∞, mle_tests_α[i]= fit_mle_α(wienerARD∞, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD∞.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD∞.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD∞.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α
end

function tests_β(wienerARD∞::WienerARD∞, degradationdata::DegradationData, nb_tests::Int64, β::Float64)#version où il suffit de rentrer une structure de degradationdata et il fera les simulations; approprie pour les structures assymetriques
    mle_tests_drift = Vector{Float64}(undef, nb_tests)
    mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
    mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD∞.maintenances), nb_tests)
    mle_tests_β = Vector{Float64}(undef, nb_tests)

    for i in 1:nb_tests
        generateur_laurent_non_lin(wienerARD∞, degradationdata, 0, x -> x ^ α, x -> x ^ β)
        
        ith_wienerARD∞, mle_tests_β[i] = fit_mle_β(wienerARD∞, degradationdata)

        mle_tests_maintenance_effect[:, i] = ith_wienerARD∞.maintenances.VALUE
        mle_tests_drift[i] = ith_wienerARD∞.underlyingDegradation[1]
        mle_tests_dispersion[i] = ith_wienerARD∞.underlyingDegradation[2]^2
    end

    return mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_β
end


#####PLOT FUNCTIONS FOR SIMULATIONS#####
function plot_tests(wienerARD1::WienerARD1, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion; savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD1.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD1.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    # ylims!(-.5, 4.5)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    # ylims!(0., 12.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)
    
    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\"*savename*".png", fig)
    end

    display(fig)
end

function plot_tests_α(wienerARD1::WienerARD1, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α; α, savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD1.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD1.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    ylims!(-.5, 8.)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    ylims!(0., 20.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)

    ax_α = Axis(fig[1, 3 + nb_ρ], title = L"\alpha", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_α, zeros(nb_tests), mle_tests_α, color=:magenta)
    hlines!(α, color=:blue, linestyle=:dash)
    
    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\tests_alpha\\"*savename*".png", fig)
    end

    display(fig)
end

function plot_tests_β(wienerARD1::WienerARD1, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_β; β, savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD1.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD1.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    ylims!(-.5, 8.)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    ylims!(0., 20.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)
    
    ax_β = Axis(fig[1, 3 + nb_ρ], title = L"\beta", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_β, zeros(nb_tests), mle_tests_β, color=:magenta)
    hlines!(β, color=:blue, linestyle=:dash)

    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\tests_alpha_beta\\"*savename*".png", fig)
    end

    display(fig)
end

function plot_tests_α_β(wienerARD1::WienerARD1, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α, mle_tests_β; α, β, savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD1.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD1.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    ylims!(-.5, 8.)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    ylims!(0., 20.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD1.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)

    ax_α = Axis(fig[1, 3 + nb_ρ], title = L"\alpha", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_α, zeros(nb_tests), mle_tests_α, color=:magenta)
    hlines!(α, color=:blue, linestyle=:dash)
    
    ax_β = Axis(fig[1, 4 + nb_ρ], title = L"\beta", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_β, zeros(nb_tests), mle_tests_β, color=:magenta)
    hlines!(β, color=:blue, linestyle=:dash)

    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\tests_alpha_beta\\"*savename*".png", fig)
    end

    display(fig)
end

function plot_tests_α(wienerARD∞::WienerARD∞, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α; α, savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD∞.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD∞.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    ylims!(-.5, 8.)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    ylims!(0., 20.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)

    ax_α = Axis(fig[1, 3 + nb_ρ], title = L"\alpha", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_α, zeros(nb_tests), mle_tests_α, color=:magenta)
    hlines!(α, color=:blue, linestyle=:dash)
    
    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\tests_alpha\\"*savename*".png", fig)
    end

    display(fig)
end

function plot_tests_β(wienerARD∞::WienerARD∞, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α, mle_tests_β; β, savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD∞.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD∞.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    ylims!(-.5, 8.)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    ylims!(0., 20.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)
    
    ax_β = Axis(fig[1, 3 + nb_ρ], title = L"\beta", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_β, zeros(nb_tests), mle_tests_β, color=:magenta)
    hlines!(β, color=:blue, linestyle=:dash)

    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\tests_alpha_beta\\"*savename*".png", fig)
    end

    display(fig)
end

function plot_tests_α_β(wienerARD∞::WienerARD∞, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion, mle_tests_α, mle_tests_β; α, β, savepicture=false, savename="picture")
    fig = Figure()

    nb_ρ = nrow(wienerARD∞.maintenances)
    for i in 1:nb_ρ
        ax_maintenance_effect = Axis(fig[1, i], title = L"\rho_{%$i}", titlesize=20., xticks = ([], []))
        ylims!(-.25, 1.25)
        boxplot!(ax_maintenance_effect, fill(i, nb_tests), mle_tests_maintenance_effect[i, :], color=:magenta)
        hlines!(wienerARD∞.maintenances[i, "VALUE"], color=:blue, linestyle=:dash)
    end

    ax_drift = Axis(fig[1, 1 + nb_ρ], title = L"\mu", titlesize=20, xticks = ([], []))
    ylims!(-.5, 8.)
    boxplot!(ax_drift, ones(nb_tests), mle_tests_drift, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[1], color=:blue, linestyle=:dash)
    
    ax_dispersion = Axis(fig[1, 2 + nb_ρ], title = L"\sigma^2", titlesize=20, xticks = ([], []))
    ylims!(0., 20.)
    boxplot!(ax_dispersion, zeros(nb_tests), mle_tests_dispersion, color=:magenta)
    hlines!(wienerARD∞.underlyingDegradation[2]^2, color=:blue, linestyle=:dash)

    ax_α = Axis(fig[1, 3 + nb_ρ], title = L"\alpha", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_α, zeros(nb_tests), mle_tests_α, color=:magenta)
    hlines!(α, color=:blue, linestyle=:dash)
    
    ax_β = Axis(fig[1, 4 + nb_ρ], title = L"\beta", titlesize=20, xticks = ([], []))
    ylims!(0., 2.)
    boxplot!(ax_β, zeros(nb_tests), mle_tests_β, color=:magenta)
    hlines!(β, color=:blue, linestyle=:dash)

    if savepicture == true
        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\pictures\\tests_alpha_beta\\"*savename*".png", fig)
    end

    display(fig)
end