using Pkg
Pkg.activate(".")

using Revise
using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

include("unidim_functests_and_plottests.jl")

using CSV
using DataFrames
using Random

#####EDF DATA 1#####
    deg1 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\deg1.csv", DataFrame)
    maint1 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\maint1.csv", DataFrame)
    deg1.TYPE = fill("BETWEEN", nrow(deg1))
    maint1.TYPE = convert.(String, maint1.TYPE)
    degradationdata1 = DegradationData(maint1, deg1)
    deg2 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\deg2.csv", DataFrame)
    maint2 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\maint2.csv", DataFrame)
    deg2.TYPE = fill("BETWEEN", nrow(deg2))
    maint2.TYPE = convert.(String, maint2.TYPE)
    degradationdata2 = DegradationData(maint2, deg2)
    deg3 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\deg3.csv", DataFrame)
    maint3 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\maint3.csv", DataFrame)
    deg3.TYPE = fill("BETWEEN", nrow(deg3))
    maint3.TYPE = convert.(String, maint3.TYPE)
    degradationdata3 = DegradationData(maint3, deg3)
    deg4 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\deg4.csv", DataFrame)
    maint4 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data1\\maint4.csv", DataFrame)
    deg4.TYPE = fill("BETWEEN", nrow(deg4))
    maint4.TYPE = convert.(String, maint4.TYPE)
    degradationdata4 = DegradationData(maint4, deg4)
    multiple_degradationdata = [degradationdata1, degradationdata2, degradationdata3, degradationdata4]
    wienerARD1 = WienerARD1(["p"], 0., 0., [.5])
    wienerARD∞ = WienerARD∞(["p"], 0., 0., [.5])

fit_mle(wienerARD1, multiple_degradationdata)
fit_mle(wienerARD1, degradationdata1)
fit_mle(wienerARD1, degradationdata2)
fit_mle(wienerARD1, degradationdata3)
fit_mle(wienerARD1, degradationdata4)

fit_mle(wienerARD∞, multiple_degradationdata)
fit_mle(wienerARD∞, degradationdata1)
fit_mle(wienerARD∞, degradationdata2)
fit_mle(wienerARD∞, degradationdata3)
fit_mle(wienerARD∞, degradationdata4)

fit_mle_α(wienerARD1, degradationdata1)
fit_mle_α(wienerARD1, degradationdata2)
fit_mle_α(wienerARD1, degradationdata3)
fit_mle_α(wienerARD1, degradationdata4)

fit_mle_α_β(wienerARD1, degradationdata1)
fit_mle_α_β(wienerARD1, degradationdata2)
fit_mle_α_β(wienerARD1, degradationdata3)
fit_mle_α_β(wienerARD1, degradationdata4)

fit_mle_α(wienerARD∞, degradationdata1)
fit_mle_α(wienerARD∞, degradationdata2)
fit_mle_α(wienerARD∞, degradationdata3)
fit_mle_α(wienerARD∞, degradationdata4)
fit_mle_α(wienerARD∞, multiple_degradationdata)

fit_mle_α_β(wienerARD∞, degradationdata1)
fit_mle_α_β(wienerARD∞, degradationdata2)
fit_mle_α_β(wienerARD∞, degradationdata3)
fit_mle_α_β(wienerARD∞, degradationdata4)
fit_mle_α_β(wienerARD∞, multiple_degradationdata)

multiple_Y, multiple_E, multiple_V = coefficient([.5], multiple_degradationdata, wienerARD∞, 1., 1.)
drift = fit_mle_drift∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, 1., 1.)
fit_mle_dispersion∞(multiple_degradationdata, multiple_Y, multiple_E, multiple_V, drift, 1., 1.)
fit_mle_maintenance_effect_α(wienerARD∞, multiple_degradationdata).minimizer[1:2]
fit_mle_α(wienerARD∞, multiple_degradationdata)
fit_mle_maintenance_effect_α_β(wienerARD∞, multiple_degradationdata).minimizer
fit_mle_α_β(wienerARD∞, multiple_degradationdata)

#####FIGURES TO PRESENT THE TRAJECTORIES#####
    fig = Figure()
    ax1, ax2, ax3, ax4 = Axis(fig[1, 1]), Axis(fig[1, 2]), Axis(fig[2, 1]), Axis(fig[2, 2])
    plot!(ax1, degradationdata1)
    plot!(ax2, degradationdata2)
    plot!(ax3, degradationdata3)
    plot!(ax4, degradationdata4)
    display(fig)

    fig = Figure()
    ax = Axis(fig[1, 1], xticks = ([], []))
    plot!(ax, degradationdata3)
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\EDFdata1_3.png", fig)

    fig = Figure()
    ax = Axis(fig[1, 1], xticks = ([], []))
    plot!(ax, degradationdata4)
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\EDFdata1_4.png", fig)

nb_tests = 5000
mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion = tests(wienerARD1, multiple_degradationdata, nb_tests)
plot_tests(wienerARD1, mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion)



#####EDF DATA 2#####
    ###NAMING###
    deg1 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg1.csv", DataFrame)
    maint1 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint1.csv", DataFrame)
    deg1.TYPE = fill("BETWEEN", nrow(deg1))
    maint1.TYPE = convert.(String, maint1.TYPE)
    degradationdata1 = DegradationData(maint1, deg1)
    deg2 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg2.csv", DataFrame)
    maint2 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint2.csv", DataFrame)
    deg2.TYPE = fill("BETWEEN", nrow(deg2))
    maint2.TYPE = convert.(String, maint2.TYPE)
    degradationdata2 = DegradationData(maint2, deg2)
    deg3 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg3.csv", DataFrame)
    maint3 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint3.csv", DataFrame)
    deg3.TYPE = fill("BETWEEN", nrow(deg3))
    maint3.TYPE = convert.(String, maint3.TYPE)
    degradationdata3 = DegradationData(maint3, deg3)
    deg4 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg4.csv", DataFrame)
    maint4 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint4.csv", DataFrame)
    deg4.TYPE = fill("BETWEEN", nrow(deg4))
    maint4.TYPE = convert.(String, maint4.TYPE)
    degradationdata4 = DegradationData(maint4, deg4)
    deg5 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg5.csv", DataFrame)
    maint5 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint5.csv", DataFrame)
    deg5.TYPE = fill("BETWEEN", nrow(deg5))
    maint5.TYPE = convert.(String, maint5.TYPE)
    degradationdata5 = DegradationData(maint5, deg5)
    deg6 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg6.csv", DataFrame)
    maint6 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint6.csv", DataFrame)
    deg6.TYPE = fill("BETWEEN", nrow(deg6))
    maint6.TYPE = convert.(String, maint6.TYPE)
    degradationdata6 = DegradationData(maint6, deg6)
    deg7 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg7.csv", DataFrame)
    maint7 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint7.csv", DataFrame)
    deg7.TYPE = fill("BETWEEN", nrow(deg7))
    maint7.TYPE = convert.(String, maint7.TYPE)
    degradationdata7 = DegradationData(maint7, deg7)
    deg8 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg8.csv", DataFrame)
    maint8 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint8.csv", DataFrame)
    deg8.TYPE = fill("BETWEEN", nrow(deg8))
    maint8.TYPE = convert.(String, maint8.TYPE)
    degradationdata8 = DegradationData(maint8, deg8)
    deg9 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\deg9.csv", DataFrame)
    maint9 = CSV.read("C:\\Users\\bergerto\\Documents\\PhD\\Julia\\ProgJuliav2\\ARDmodels\\exemples\\EDF_data2\\maint9.csv", DataFrame)
    deg9.TYPE = fill("BETWEEN", nrow(deg9))
    maint9.TYPE = convert.(String, maint9.TYPE)
    degradationdata9 = DegradationData(maint9, deg9)
    multiple_degradationdata = [degradationdata1, degradationdata2, degradationdata3, degradationdata4, degradationdata5, degradationdata6, degradationdata7, degradationdata8, degradationdata9]
    maintTypes = [ "P0", "P1", "P2"]
    wienerARD1 = WienerARD1(maintTypes, 0., 0., [.5, .5, .5])
    wienerARD∞ = WienerARD∞(maintTypes, 0., 0., [.5, .5, .5])

    #####ESTIMATION#####
        ###ARD1###
            fit_mle(wienerARD1, multiple_degradationdata)
            fit_mle(wienerARD1, degradationdata1)
            fit_mle(wienerARD1, degradationdata2)
            fit_mle(wienerARD1, degradationdata3)
            fit_mle(wienerARD1, degradationdata4)
            fit_mle(wienerARD1, degradationdata5)
            fit_mle(wienerARD1, degradationdata6)
            fit_mle(wienerARD1, degradationdata7)
            fit_mle(wienerARD1, degradationdata8)
            fit_mle(wienerARD1, degradationdata9)

            ARD.fit_mle_α_β(wienerARD1, multiple_degradationdata)
            ARD.fit_mle_α(wienerARD1, multiple_degradationdata)
            ARD.fit_mle_β(wienerARD1, multiple_degradationdata)


        ###ARD∞###
            fit_mle(wienerARD∞, multiple_degradationdata)
            fit_mle(wienerARD∞, degradationdata1)
            fit_mle(wienerARD∞, degradationdata2)
            fit_mle(wienerARD∞, degradationdata3)
            fit_mle(wienerARD∞, degradationdata4)
            fit_mle(wienerARD∞, degradationdata5)
            fit_mle(wienerARD∞, degradationdata6)
            fit_mle(wienerARD∞, degradationdata7)
            fit_mle(wienerARD∞, degradationdata8)
            fit_mle(wienerARD∞, degradationdata9)

            ARD.fit_mle_α_β(wienerARD∞, multiple_degradationdata)
            ARD.fit_mle_α(wienerARD∞, multiple_degradationdata)
            ARD.fit_mle_β(wienerARD∞, multiple_degradationdata)

            nb_systeme = length(multiple_degradationdata)
            multiple_Y = Vector{Vector{Float64}}(undef, nb_systeme)
            multiple_E = Vector{Vector{Float64}}(undef, nb_systeme)
            multiple_V = Vector{Vector{Float64}}(undef, nb_systeme)
            
            for i in eachindex(multiple_degradationdata)
                multiple_Y[i], multiple_E[i], multiple_V[i] = coefficient(ρ, multiple_degradationdata[1], wienerARD∞, α, β)[1:3]
            end



#####FIGURES TO PRESENT THE TRAJECTORIES#####
    for (index, degradationdata) in enumerate(multiple_degradationdata)
        fig = Figure()
        ax1 = Axis(fig[1, 1])
        plot!(ax1, degradationdata, maintenanceslinestyles=DataFrame(TYPE = maintTypes, LINESTYLE = [:solid, :dash, :dot]))
        xlims!(0., 240000.)
        ylims!(0., 35.)

        hidexdecorations!(ax1, ticklabels = true, label = false)
        hideydecorations!(ax1, ticklabels = true, label = false)

        save("C:\\Users\\bergerto\\Documents\\Stage_M2_recherche\\Julia\\Tom_Berger\\mle_gap_maintenance\\EDFdata2_$index.png", fig)
    end

#####SIMULATIONS USING THE COMBINATION OF THE 1ST, 4TH AND 5TH EDFdata2 SCHEMES OF OBSERVATION#####
using StatsBase
using Random

Random.seed!(1)
nb_tests = 5000
mle_tests_maintenance_effect, mle_tests_drift, mle_tests_dispersion = tests(wienerARD1, multiple_degradationdata, nb_tests)
mle_tests_maintenance_effect2, mle_tests_drift2, mle_tests_dispersion2 = tests(wienerARD1, vcat(multiple_degradationdata, multiple_degradationdata), nb_tests)

mse_simulation_application_df = DataFrame(":drift" => [msd(mle_tests_drift, [wienerARD1.underlyingDegradation[1] for _ in 1:nb_tests]), msd(mle_tests_drift2, [wienerARD1.underlyingDegradation[1] for _ in 1:nb_tests])],
    ":volatility" => [msd(mle_tests_dispersion, [wienerARD1.underlyingDegradation[2] for _ in 1:nb_tests]), msd(mle_tests_dispersion2, [wienerARD1.underlyingDegradation[2] for _ in 1:nb_tests])],
    ":maint1" => [msd(mle_tests_maintenance_effect[1, :], [wienerARD1.maintenances.VALUE[1] for _ in 1:nb_tests]), msd(mle_tests_maintenance_effect2[1, :], [wienerARD1.maintenances.VALUE[1] for _ in 1:nb_tests])],
    ":maint2" => [msd(mle_tests_maintenance_effect[2, :], [wienerARD1.maintenances.VALUE[2] for _ in 1:nb_tests]), msd(mle_tests_maintenance_effect2[2, :], [wienerARD1.maintenances.VALUE[2] for _ in 1:nb_tests])],
    ":maint3" => [msd(mle_tests_maintenance_effect[3, :], [wienerARD1.maintenances.VALUE[3] for _ in 1:nb_tests]), msd(mle_tests_maintenance_effect2[3, :], [wienerARD1.maintenances.VALUE[3] for _ in 1:nb_tests])])

pretty_table(mse_simulation_application_df; backend=:latex)



#####SIMULATIONS USING THE COMBINATION OF ALL EDFdata2 SCHEMES OF OBSERVATION#####
A1 = ARD.fit_mle_α(wienerARD1, multiple_degradationdata)
A∞ = ARD.fit_mle_α(wienerARD∞, multiple_degradationdata)
nlwienerard1 = ARD.NLWienerARD1(A1[1].underlyingDegradation[1], A1[1].underlyingDegradation[2], A1[2], 1., Dict("P0" => A1[1].maintenances[1, 2], "P1" => A1[1].maintenances[2, 2], "P2" => A1[1].maintenances[3, 2]))
nlwienerard∞ = ARD.NLWienerARD∞(A∞[1].underlyingDegradation[1], A∞[1].underlyingDegradation[2], A∞[2], 1., Dict("P0" => A∞[1].maintenances[1, 2], "P1" => A∞[1].maintenances[2, 2], "P2" => A∞[1].maintenances[3, 2]))

Random.seed!(1)
nb_tests = 5000

mle_tests_drift = Vector{Float64}(undef, nb_tests)
mle_tests_dispersion = Vector{Float64}(undef, nb_tests)
mle_tests_maintenance_effect = Array{Float64}(undef, nrow(wienerARD1.maintenances), nb_tests)
mle_tests_α = Vector{Float64}(undef, nb_tests)

###ARD1 case###
for i in 1:nb_tests
    for degradationdata in multiple_degradationdata
        generateur_laurent_non_lin(A1[1], degradationdata, 0, x -> x ^ A1[2], x -> x)
    end
    
    ith_wienerARD1, mle_tests_α[i]= ARD.fit_mle_α(wienerARD1, multiple_degradationdata)

    mle_tests_maintenance_effect[:, i] = ith_wienerARD1.maintenances.VALUE
    mle_tests_drift[i] = ith_wienerARD1.underlyingDegradation[1]
    mle_tests_dispersion[i] = ith_wienerARD1.underlyingDegradation[2]^2
end

using StatsBase

mse_simulation_application_df = DataFrame(":drift" => [msd(mle_tests_drift, [wienerARD1.underlyingDegradation[1] for _ in 1:nb_tests])],
    ":volatility" => [msd(mle_tests_dispersion, [wienerARD1.underlyingDegradation[2] for _ in 1:nb_tests])],
    ":maint1" => [msd(mle_tests_maintenance_effect[1, :], [wienerARD1.maintenances.VALUE[1] for _ in 1:nb_tests])],
    ":maint2" => [msd(mle_tests_maintenance_effect[2, :], [wienerARD1.maintenances.VALUE[2] for _ in 1:nb_tests])],
    ":maint3" => [msd(mle_tests_maintenance_effect[3, :], [wienerARD1.maintenances.VALUE[3] for _ in 1:nb_tests])],
    ":α" => [msd(mle_tests_α, [A1[2] for _ in 1:nb_tests])])

using PrettyTables

pretty_table(mse_simulation_application_df; backend=:latex)

###ARD∞ case###
for i in 1:nb_tests
    for degradationdata in multiple_degradationdata
        generateur_laurent_non_lin(A∞[1], degradationdata, 0, x -> x ^ A∞[2], x -> x)
    end
    
    ith_wienerARD∞, mle_tests_α[i]= ARD.fit_mle_α(wienerARD∞, multiple_degradationdata)

    mle_tests_maintenance_effect[:, i] = ith_wienerARD∞.maintenances.VALUE
    mle_tests_drift[i] = ith_wienerARD∞.underlyingDegradation[1]
    mle_tests_dispersion[i] = ith_wienerARD∞.underlyingDegradation[2]^2
end

using StatsBase

mse_simulation_application_df = DataFrame(":drift" => [msd(mle_tests_drift, [wienerARD∞.underlyingDegradation[1] for _ in 1:nb_tests])],
    ":volatility" => [msd(mle_tests_dispersion, [wienerARD∞.underlyingDegradation[2] for _ in 1:nb_tests])],
    ":maint1" => [msd(mle_tests_maintenance_effect[1, :], [wienerARD∞.maintenances.VALUE[1] for _ in 1:nb_tests])],
    ":maint2" => [msd(mle_tests_maintenance_effect[2, :], [wienerARD∞.maintenances.VALUE[2] for _ in 1:nb_tests])],
    ":maint3" => [msd(mle_tests_maintenance_effect[3, :], [wienerARD∞.maintenances.VALUE[3] for _ in 1:nb_tests])],
    ":α" => [msd(mle_tests_α, [A∞[2] for _ in 1:nb_tests])])

using PrettyTables

pretty_table(mse_simulation_application_df; backend=:latex)

ARD.fit_mle_α(wienerARD1, multiple_degradationdata)
for _ in 1:nb_tests
    mle_tests_maintenance_effect[i], mle_tests_drift[i], mle_tests_dispersion[i], mle_tests_α[i] = fit_mle_α(wienerARD1, multiple_degradationdata)
end