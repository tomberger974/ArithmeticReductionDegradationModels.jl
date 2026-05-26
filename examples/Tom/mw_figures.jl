using Latexify
using CairoMakie
using DataFrames
using Distributions
using Random
using LinearAlgebra

using Pkg
Pkg.activate(".")

using Revise
using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

include("nearest_time.jl")

# Define the multilinear ARD processes
mward1 = ARD.MWARD1([2., 2.], [.4 .3; .3 .4], Dict("P" => [0.5, 0.2], "C" => [0.8, 0.5], "T" => [.9, .9]))
mward∞ = ARD.MWARD∞([2., 2.], [.4 .3; .3 .4], Dict("P" => [0.5, 0.2], "C" => [0.8, 0.5], "T" => [.9, .9]))
[.3 .5; .5 .3]
# Define the inspection dates
T = 15.
inspection_dates = convert(Vector{Float64}, range(0, T, 500))

# Maintenance dates
τ = T .* [.17, .4, .5, .8]
maintenances = DataFrame(DATE=τ, TYPE=["P", "C", "P", "T"])

#if there is only one time stop τ there is no difference between ARD1 and ARD∞
Random.seed!(3)
Y1 = rand(mward1, inspection_dates, maintenances)
Random.seed!(3)
Y∞ = rand(mward∞, inspection_dates, maintenances)




##### TWO CORRELATED UNMAINTAINED PROCESS ######
# Inspection dates indices
I = [18, 56, 102, 270, 410, 464, 498]
I_char1 = I[[3, 4, 5, 7]]
I_char2 = I[[1, 2, 6, 7]]

# LINEAR WIENER PROCESS
fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1], ylabel = "(1)", xgridvisible=false)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Time", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1], ylabel = "(2)", xgridvisible=false)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:cyan)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple)
    vlines!(ax1, τ[3], linestyle=:dash, color=:cyan)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:orange)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:cyan)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple)
    vlines!(ax2, τ[3], linestyle=:dash, color=:cyan)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:orange)

    # Main Wiener process trajectory
    Random.seed!(3)
    main_trajectory = ARD.rand(mward1, inspection_dates, maintenances)
    char1 = lines!(ax1, inspection_dates, main_trajectory[1, :], color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, main_trajectory[2, :], color=:blue, linewidth=1.)

    # Drift directing curve
    drift_char1 = lines!(ax1, inspection_dates, mward1.drift[1] .* inspection_dates, color=:blue, alpha=0.5)
    drift_char2 = lines!(ax2, inspection_dates, mward1.drift[2] .* inspection_dates, color=:blue, alpha=0.5)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[1, I_char1], color=:red, markersize=15.)
    vlines!(ax1, inspection_dates[I_char1], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y1[2, I_char2], color=:red, markersize=15.)
    vlines!(ax2, inspection_dates[I_char2], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=20.)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=20.)

    # Display and save the figure
    display(fig)




##### FIVE MAINTENANCE ACTIONS #####
# NON-LINEAR WIENER PROCESS WITH SQUARE ROOT DRIFT
fig = Figure(resolution = (1920, 1080))

    # Axis definition
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])
    ax3 = Axis(fig[3, 1])

    # Make the second figure smaller
    rowsize!(fig.layout, 2, Relative(0.1))

    # Main Wiener process trajectory
    Random.seed!(7)
    main_trajectory = ARD.rand(mward1, inspection_dates, maintenances)
    char1 = lines!(ax1, inspection_dates, main_trajectory[1, :], color=:blue, linewidth=1.)
    char2 = lines!(ax3, inspection_dates, main_trajectory[2, :], color=:red, linewidth=1.)

    # Secondary trajectories to illustrate the variability of the process
    for i in 1:5
        traj = ARD.mw_rand(mward1, inspection_dates)
        lines!(ax1, inspection_dates, traj[1, :], color = (:grey, 0.5), linewidth=1.)
        lines!(ax3, inspection_dates, traj[2, :], color = (:grey, 0.5), linewidth=1.)
    end

    # Drift directing curve
    # droite = lines!(ax1, inspection_dates, nlwienerard1.drift .* inspection_dates .^ nlwienerard1.power_drift)

    # Display and save the figure
    display(fig)