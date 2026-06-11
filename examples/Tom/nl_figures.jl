using Latexify
using CairoMakie
using DataFrames
using Distributions
using Random

using Pkg
Pkg.activate(".")

using Revise
using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

include("nearest_time.jl")

nlwienerard1 = ARD.NLWienerARD1(2., .25, .5, 1., Dict("P" => .4, "C" => .9, "T" => .6))
nlwienerard∞ = ARD.NLWienerARD∞(2., .25, .5, 1., Dict("P" => .4, "C" => .9, "T" => .6))


# Inspection dates
T = 15.
inspection_dates = convert(Vector{Float64}, range(0, T, 500))

# Maintenance dates
τ = T .* [.17, .4, .5, .7, .9]
maintenances = DataFrame(DATE=τ, TYPE=["P", "C", "P", "P", "T"])

#if there is only one time stop τ there is no difference between ARD1 and ARD∞
Random.seed!(7)
Y1 = rand(nlwienerard∞, inspection_dates, maintenances)
Random.seed!(7)
Y∞ = rand(nlwienerard1, inspection_dates, maintenances)




#### NO MAINTENANCE #####


# NON-LINEAR WIENER PROCESS WITH SQUARE ROOT DRIFT
fig = Figure()

    # Axis definition
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Degradation Value")

    # Main Wiener process trajectory
    Random.seed!(7)
    wiener_process_plot = lines!(ax, inspection_dates, ARD.nl_wiener_rand(nlwienerard1, inspection_dates), color=:blue, linewidth=1.)

    # Secondary trajectories to illustrate the variability of the process
    for i in 1:5
        traj = ARD.nl_wiener_rand(nlwienerard1, inspection_dates)
        lines!(ax, inspection_dates, traj, color = (:grey, 0.5), linewidth=1.)
    end

    # Drift directing curve
    droite = lines!(ax, inspection_dates, nlwienerard1.drift .* inspection_dates .^ nlwienerard1.power_drift)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_wiener_process_multiple_traj_sqrt.png", fig)


# LINEAR WIENER PROCESS
fig = Figure()

    # Axis definition
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Degradation Value")

    # Main Wiener process trajectory
    Random.seed!(7)
    wiener_process_plot = lines!(ax, inspection_dates, ARD.nl_wiener_rand(ARD.NLWienerARD1(2., 1., 1., 1., Dict("P" => .4, "C" => .9, "T" => .6)), inspection_dates), color=:blue, linewidth=1.)

    # Secondary trajectories to illustrate the variability of the process
    for i in 1:5
        traj = ARD.nl_wiener_rand(ARD.NLWienerARD1(2., 1., 1., 1., Dict("P" => .4, "C" => .9, "T" => .6)), inspection_dates)
        lines!(ax, inspection_dates, traj, color = (:grey, 0.5), linewidth=1.)
    end

    # Drift directing curve
    droite = lines!(ax, inspection_dates, nlwienerard1.drift .* inspection_dates)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_wiener_process_multiple_traj_linear.png", fig)




#### FIVE MAINTENANCE, NON-LINEAR WIENER PROCESS WITH SQUARE ROOT DRIFT #####


# BOTH ARD1 AND ARD∞, ILLUSTRATION OF ARD FUNCTIONMENT
fig = Figure()

    # Axis definition with an axis for the inspection dates and another one for the maintenance dates
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "Time", ylabel = "Degradation Value")
    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(τ[1], linestyle=:dash, color=:cyan)
    maint_type2 = vlines!(τ[2], linestyle=:dash, color=:purple)
    vlines!(τ[3], linestyle=:dash, color=:cyan)
    vlines!(τ[4], linestyle=:dash, color=:cyan)
    maint_type3 = vlines!(τ[5], linestyle=:dash, color=:orange)

    # Plot ARD trajectories
    ARD∞_plot = lines!(ax, inspection_dates, Y1, color=:red, linewidth=1.)
    ARD1_plot = lines!(ax, inspection_dates, Y∞, color=:blue, linewidth=1.)

    # ARD1 arrows
    arrows!(ax, (inspection_dates[nearest_time(τ[2], inspection_dates) - 10], Y∞[nearest_time(τ[2], inspection_dates) - 1]), (0., Y∞[nearest_time(τ[2], inspection_dates)] - Y∞[nearest_time(τ[2], inspection_dates) - 1]), color = :blue)
    text!(inspection_dates[nearest_time(τ[2], inspection_dates) - 55], (Y∞[nearest_time(τ[2], inspection_dates) - 1] + Y∞[nearest_time(τ[2], inspection_dates)]) / 2 - .5, text=L"\rho_2 x", color=:blue, fontsize=20.)
    arrows!(ax, (inspection_dates[nearest_time(τ[1], inspection_dates) + 20], Y∞[nearest_time(τ[1], inspection_dates)]), (0., Y∞[nearest_time(τ[2], inspection_dates) - 1] - Y∞[nearest_time(τ[1], inspection_dates)]), color = :blue)
    text!(inspection_dates[nearest_time(τ[1], inspection_dates) + 20], (Y∞[nearest_time(τ[1], inspection_dates)] + Y∞[nearest_time(τ[2], inspection_dates) - 1]) / 2, text=L"x", color=:blue, fontsize=20.)

    # ARD∞ arrows
    arrows!(ax, (inspection_dates[nearest_time(τ[2], inspection_dates) - 20], Y1[nearest_time(τ[2], inspection_dates) - 1]), (0., - nlwienerard1.efficiencies["C"] * Y1[nearest_time(τ[2], inspection_dates) - 1]), color = :red)
    text!(inspection_dates[nearest_time(τ[2], inspection_dates) - 55], (Y1[nearest_time(τ[2], inspection_dates) - 1] + Y1[nearest_time(τ[2], inspection_dates)]) / 2 - 1, text=L"\rho_2 x", color=:red, fontsize=20.)
    arrows!(ax, (inspection_dates[nearest_time(τ[1], inspection_dates) + 10], 0.), (0., Y1[nearest_time(τ[2], inspection_dates) - 1]), color = :red)
    text!(inspection_dates[nearest_time(τ[1], inspection_dates) + 10], Y1[nearest_time(τ[2], inspection_dates) - 1] / 4, text=L"x", color=:red, fontsize=20.)

    # Grey dashed lines to illustrate the values of the process at maintenance dates
    lines!(inspection_dates[[1, nearest_time(τ[2], inspection_dates)]], [Y∞[nearest_time(τ[2], inspection_dates) - 1] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(inspection_dates[[nearest_time(τ[1], inspection_dates), nearest_time(τ[2], inspection_dates)]], [Y1[nearest_time(τ[1], inspection_dates)] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(inspection_dates[[1, nearest_time(τ[2], inspection_dates)]], [0. for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)

    # Legend
    axislegend(ax, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1 = 0.4", L"\rho_2 = 0.9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_both_degradation_process.png", fig)


# ARD1 WITH ILLUSTRATION OF ARD FUNCTIONMENT ON THE FIFTH MAINTENANCE ACTION
fig = Figure()

    # Axis definition with an axis for the inspection dates and another one for the maintenance dates
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "Time", ylabel = "Degradation Value")
    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(τ[1], linestyle=:dash, color=:cyan)
    maint_type2 = vlines!(τ[2], linestyle=:dash, color=:purple)
    vlines!(τ[3], linestyle=:dash, color=:cyan)
    vlines!(τ[4], linestyle=:dash, color=:cyan)
    maint_type3 = vlines!(τ[5], linestyle=:dash, color=:orange)

    # Plot ARD trajectories
    ARD1_plot = lines!(ax, inspection_dates, Y∞, color=:blue, linewidth=1.)

    # Arrows to illustrate the effect of the 5th maintenance on ARD1
    arrows!(ax, (inspection_dates[nearest_time(τ[5], inspection_dates) - 10], Y∞[nearest_time(τ[5], inspection_dates) - 1]), (0., Y∞[nearest_time(τ[5], inspection_dates)] - Y∞[nearest_time(τ[5], inspection_dates) - 1]), color = :orange)
    text!(inspection_dates[nearest_time(τ[5], inspection_dates) - 60], (Y∞[nearest_time(τ[5], inspection_dates)] + Y∞[nearest_time(τ[5], inspection_dates) - 1]) / 2, text=L"\rho_3 x", color=:orange, fontsize=20.)
    arrows!(ax, (inspection_dates[nearest_time(τ[5], inspection_dates) - 60], Y∞[nearest_time(τ[4], inspection_dates)]), (0., Y∞[nearest_time(τ[5], inspection_dates) - 1] - Y∞[nearest_time(τ[4], inspection_dates)]), color = :purple)
    text!(inspection_dates[nearest_time(τ[5], inspection_dates) - 75], (Y∞[nearest_time(τ[4], inspection_dates)] + Y∞[nearest_time(τ[5], inspection_dates) - 1]) / 2, text=L"x", color=:purple, fontsize=20.)

    # Grey dashed lines to illustrate the values of the process at maintenance dates
    lines!(inspection_dates[[nearest_time(τ[4], inspection_dates), nearest_time(τ[5], inspection_dates)]], [Y∞[nearest_time(τ[5], inspection_dates) - 1] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(inspection_dates[[nearest_time(τ[4], inspection_dates), nearest_time(τ[5], inspection_dates)]], [Y∞[nearest_time(τ[5], inspection_dates)] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(inspection_dates[[nearest_time(τ[4], inspection_dates), nearest_time(τ[5], inspection_dates)]], [Y∞[nearest_time(τ[4], inspection_dates)] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)

    # Legend
    #axislegend(ax, merge = true, [ARD1_plot, ARD∞_plot, maint_type1, maint_type2, maint_type3], [L"ARD_1", L"ARD_\infty", L"\rho_1 = 0.4", L"\rho_2 = .9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)

    # Display the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_degradation_process_how_ard1_works.png", fig)


# ILLUSTRATES SAMPLING, ADD OF THE OBSERVATION POINTS, ONLY ARD1
fig = Figure()

    # Axis definition with an axis for the inspection dates and another one for the maintenance dates
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "", ylabel = "Degradation Value")
    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    # Maintenance datesillustration with dashed colored vertical lines
    maint_type1 = vlines!(τ[1], linestyle=:dash, color=:cyan)
    maint_type2 = vlines!(τ[2], linestyle=:dash, color=:purple)
    vlines!(τ[3], linestyle=:dash, color=:cyan)
    vlines!(τ[4], linestyle=:dash, color=:cyan)
    maint_type3 = vlines!(τ[5], linestyle=:dash, color=:orange)

    # Plot ARD1 trajectory
    ARD1_plot = lines!(ax, inspection_dates, Y∞, linewidth=1., color=:blue)

    # Add observation points with red color and bigger markersize
    I = [56, 102, 190, 270, 464, 498]
    observations = scatter!(ax, inspection_dates[I], Y∞[I], color=:red, markersize=15.)

    # Legend
    axislegend(ax, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1 = 0.4", L"\rho_2=0.9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_degradation_process_plus_observations.png", fig)


# ILLUSTRATES SAMPLING, ADD OF THE OBSERVATION POINTS, ONLY ARD1, WITH THE INCREMENT AND THE JUMP ILLUSTRATION
fig = Figure()

    # Axis definition with an axis for the inspection dates and another one for the maintenance dates
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "Time", ylabel = "Degradation Value")
    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(τ[1], linestyle=:dash, color=:cyan)
    maint_type2 = vlines!(τ[2], linestyle=:dash, color=:purple)
    vlines!(τ[3], linestyle=:dash, color=:cyan)
    vlines!(τ[4], linestyle=:dash, color=:cyan)
    maint_type3 = vlines!(τ[5], linestyle=:dash, color=:orange)

    # Plot ARD1 trajectory
    ARD1_plot = lines!(ax, inspection_dates, Y∞, linewidth=1., color=:blue)

    # Add observation points with red color and bigger markersize
    I = [56, 102, 190, 270, 464, 498]
    observations = scatter!(ax, inspection_dates[I], Y∞[I], color=:red, markersize=15.)

    # Usual Wiener increments
    lines!(inspection_dates[[102, 190]], Y∞[[190, 190]], linestyle=:dot, color=:grey, alpha=:.5)
    bracket!(ax, (inspection_dates[102], Y∞[102]), (inspection_dates[102], Y∞[190]), color = :red, rotation = 0)
    text!(inspection_dates[30], (Y∞[190] + Y∞[102])/2, text=L"\Delta Y_{\nu_2, 1}", color=:red, fontsize=18.)

    # Jumps
    lines!(inspection_dates[[270, 190]], Y∞[[190, 190]], linestyle=:dot, color=:grey, alpha=:.5)
    bracket!(ax, (inspection_dates[270], Y∞[190]), (inspection_dates[270], Y∞[270]), color = RGBAf(52/255, 119/255, 37/255, 1.), rotation = 0)
    text!(inspection_dates[300], (Y∞[190] + Y∞[270])/2 - 0.4, text=L"Z_{\nu_2}", color = RGBAf(52/255, 119/255, 37/255, 1.), fontsize=18.)

    # Legend
    # axislegend(ax, merge = true, [ARD1_plot, maint_type1, maint_type2, maint_type3, observations], [L"ARD_1", L"\rho_1 = 0.4", L"\rho_2=1.", L"\rho_3 = 0.6", L"Observations"], position = :lt, labelsize=20.)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\degradation_process_plus_observations_illustration_incrment_plus_jump.png", fig)




#### NO MAINTENANCE, DIFFERENT NON LINEARITIES #####

# Define a set of value for the different parameters of the non linear wiener processes
α_list = [.5, 1., 1.5]
β_list = [.5, 1., 1.5]
σ_list = [.25, 1., 4.]


# Plot the trajectories of Wiener processes with non linear drifts (and different volatilities)
fig = Figure(resolution = (1000, 400))

    for i in 1:3
        # Define the non linear wiener process with the corresponding parameters (it is an ARD process but ARD properties are not used)
        nlwienerard1 = ARD.NLWienerARD1(2., σ_list[i], α_list[i], 1., Dict("P" => .4, "C" => .9, "T" => .6))
        
        # Define the axis
        ax = Axis(fig[1, i], xlabel = "Time", ylabel = "Degradation Value")
        
        # Plot the trajectory of the non linear wiener process
        Random.seed!(7)
        wiener_process_plot = lines!(ax, inspection_dates, ARD.nl_wiener_rand(nlwienerard1, inspection_dates), color=:blue, linewidth=1.)

        # Plot some other trajectories to illustrate the variability of the process
        for _ in 1:5
            traj = ARD.nl_wiener_rand(nlwienerard1, inspection_dates)
            lines!(ax, inspection_dates, traj, color = (:grey, 0.5), linewidth=1.)
        end

        # Plot the drift directing curve
        droite = lines!(ax, inspection_dates, nlwienerard1.drift .* inspection_dates .^ nlwienerard1.power_drift)
    end

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_wiener_process_multiple_traj_multiple_alpha.png", fig)


# Plot the trajectories of Wiener processes with non linear volatility
fig = Figure(resolution = (1000, 400))
    for i in 1:3
        # Define the non linear wiener process with the corresponding parameters β (it is an ARD process but ARD properties are not used)
        nlwienerard1 = ARD.NLWienerARD1(2., 1., 1., β_list[i], Dict("P" => .4, "C" => .9, "T" => .6))
        
        # Define the axis
        ax = Axis(fig[1, i], xlabel = "Time", ylabel = "Degradation Value")
        
        # Plot the trajectory of the non linear wiener process
        Random.seed!(7)
        wiener_process_plot = lines!(ax, inspection_dates, ARD.nl_wiener_rand(nlwienerard1, inspection_dates), color=:blue, linewidth=1.)

        # Plot some other trajectories to illustrate the variability of the process
        for _ in 1:5
            traj = ARD.nl_wiener_rand(nlwienerard1, inspection_dates)
            lines!(ax, inspection_dates, traj, color = (:grey, 0.5), linewidth=1.)
        end

        # Plot the drift directing curve
        droite = lines!(ax, inspection_dates, nlwienerard1.drift .* inspection_dates .^ nlwienerard1.power_drift)
        
        # Set the same y limits for the 3 plots to better illustrate the effect of β on the non linearity of the drift
        if i!= 2
            ylims!(ax, 0, 35.)
        end
    end

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_wiener_process_multiple_traj_multiple_beta.png", fig)



# Plot the trajectories of Wiener processes with time scale change
fig = Figure(resolution = (1000, 400))
    for i in 1:3
        # Define the non linear wiener process with the corresponding parameters α = β
        nlwienerard1 = ARD.NLWienerARD1(2., 1., α_list[i], β_list[i], Dict("P" => .4, "C" => .9, "T" => .6))
        
        # Define the axis
        ax = Axis(fig[1, i], xlabel = "Time", ylabel = "Degradation Value")
        
        # Plot the trajectory of the non linear wiener process
        Random.seed!(7)
        wiener_process_plot = lines!(ax, inspection_dates, ARD.nl_wiener_rand(nlwienerard1, inspection_dates), color=:blue, linewidth=1.)
        
        # Plot some other trajectories to illustrate the variability of the process
        for _ in 1:5
            traj = ARD.nl_wiener_rand(nlwienerard1, inspection_dates)
            lines!(ax, inspection_dates, traj, color = (:grey, 0.5), linewidth=1.)
        end

        # Plot the drift directing curve
        droite = lines!(ax, inspection_dates, nlwienerard1.drift .* inspection_dates .^ nlwienerard1.power_drift)
    end

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nl_wiener_process_multiple_traj_multiple_alpha_beta.png", fig)