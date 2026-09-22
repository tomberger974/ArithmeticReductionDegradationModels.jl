using Latexify
using CairoMakie
using DataFrames

#include other functions
include("unidim_dataset_plus_simu_perso.jl")



#initialize some parameters
k = 500
T = 15

#####BROWNIAN MOTION#####
Random.seed!(5)
time, traj = Brownian_trajectory(k, T)

lines(time, traj)
display(plot)



#####WIENER PROCESS#####
Random.seed!(5)
μ = 2.
σ = sqrt(2)

wiener = Wiener(μ, σ)
time, traj = trajectory(wiener, k, T)

fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Degradation Value")
    wiener_process_plot = lines!(ax, time, traj, color=:black)
    droite = lines!(ax, time, μ .* time)
    #axislegend(ax, merge = true, [wiener_process_plot, droite], ["X(t)", L"\mu t"], position = :lt, labelsize = 20.)
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\wiener_process.png", fig)

fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Degradation Value")
    wiener_process_plot = lines!(ax, time, traj, color=:black, linewidth=1.)
    for i in 1:5
        time, traj = trajectory(wiener, k, T)
        lines!(ax, time, traj, color = (:grey, 0.5), linewidth=1.)
    end
    droite = lines!(ax, time, μ .* time)
#    axislegend(ax, merge = true, [droite], [L"\mu t"], position = :lt, labelsize=20.)
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\wiener_process_multiple_traj.png", fig)



#MAINTENANCE EFFECT, FIRST EXAMPLE:ONE MAINTENANCE, BOTH ARD1, ARD∞#####
τ = T .* [.17, .4, .5, .7, .9]
maint_types = ["P", "T", "C"]
maintenances = DataFrame(DATE=τ, TYPE=["P", "T", "P", "P", "C"])
ρ = [.4, .9, .6]

#define a ARD model with the same parameters μ, σ as before
wienerARD∞ = WienerARD∞(["P", "T", "C"], μ, σ, ρ)
wienerARD1 = WienerARD1(["P", "T", "C"], μ, σ, ρ)

#if there is only one time stop τ there is no difference between ARD1 and ARD∞
Random.seed!(5)
time, Y1 = trajectory(maintenances, wienerARD1, k, T)
Random.seed!(5)
time, Y∞ = trajectory(maintenances, wienerARD∞, k, T)



####MAINTENANCE EFFECT, FIRST EXAMPLE:ONE MAINTENANCE, ONLY ARD1#####
fig = Figure()

    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "Time", ylabel = "Degradation Value")
    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    # Plot some data on each axis
    ARD∞_plot = lines!(ax, time, Y∞, color=:red, linewidth=1.)
    ARD1_plot = lines!(ax, time, Y1, color=:blue, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[1]], linestyle=:dash, color=:magenta, linewidth=3.)
    maint_type2 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[2]], linestyle=:dash, color=:black, linewidth=3.)
    maint_type3 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[3]], linestyle=:dash, color=:green, linewidth=3.)

    # ARD∞ arrows
    arrows!(ax, (time[nearest_time(τ[2], time) - 10], Y∞[nearest_time(τ[2], time) - 1]), (0., Y∞[nearest_time(τ[2], time)] - Y∞[nearest_time(τ[2], time) - 1]), color = :red)
    text!(time[nearest_time(τ[2], time) - 55], (Y∞[nearest_time(τ[2], time) - 1] + Y∞[nearest_time(τ[2], time)]) / 2 - 3., text=L"\rho_2 x", color=:red, fontsize=20.)
    arrows!(ax, (time[nearest_time(τ[1], time) + 10], 0.), (0., Y∞[nearest_time(τ[2], time) - 1]), color = :red)
    text!(time[nearest_time(τ[1], time) + 10], Y∞[nearest_time(τ[2], time) - 1] / 4, text=L"x", color=:red, fontsize=20.)

    # ARD1 arrows
    arrows!(ax, (time[nearest_time(τ[2], time) - 20], Y1[nearest_time(τ[2], time) - 1]), (0., Y1[nearest_time(τ[2], time)] - Y1[nearest_time(τ[2], time) - 1]), color = :blue)
    text!(time[nearest_time(τ[2], time) - 55], (Y1[nearest_time(τ[2], time) - 1] + Y1[nearest_time(τ[2], time)]) / 2 - 2., text=L"\rho_2 x", color=:blue, fontsize=20.)
    arrows!(ax, (time[nearest_time(τ[1], time) + 20], Y1[nearest_time(τ[1], time)]), (0., Y1[nearest_time(τ[2], time) - 1] - Y1[nearest_time(τ[1], time)]), color = :blue)
    text!(time[nearest_time(τ[1], time) + 20], (Y1[nearest_time(τ[1], time)] + Y1[nearest_time(τ[2], time) - 1]) / 2, text=L"x", color=:blue, fontsize=20.)

    # Grey dashed lines to illustrate the values of the process at maintenance dates
    lines!(time[[1, nearest_time(τ[2], time)]], [Y∞[nearest_time(τ[2], time) - 1] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(time[[nearest_time(τ[1], time), nearest_time(τ[2], time)]], [Y1[nearest_time(τ[1], time)] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(time[[1, nearest_time(τ[2], time)]], [0. for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)

    # Legend
    axislegend(ax, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1 = 0.4", L"\rho_2 = 0.9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)
    ylims!(0., 25.)

    # Display the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\both_degradation_process.png", fig)



####MAINTENANCE EFFECT, FIRST EXAMPLE:ONE MAINTENANCE, ONLY ARD1#####
fig = Figure()

    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "Time", ylabel = "Degradation Value")
    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    # Plot some data on each axis
    ARD1_plot = lines!(ax, time, Y1, color=:blue, linewidth=1.)

    # ARD1 arrows
    arrows!(ax, (time[nearest_time(τ[2], time) - 20], Y1[nearest_time(τ[2], time) - 1]), (0., Y1[nearest_time(τ[2], time)] - Y1[nearest_time(τ[2], time) - 1]), color = :blue)
    text!(time[nearest_time(τ[2], time) - 55], (Y1[nearest_time(τ[2], time) - 1] + Y1[nearest_time(τ[2], time)]) / 2 - 2., text=L"\rho_2 x", color=:blue, fontsize=20.)
    arrows!(ax, (time[nearest_time(τ[1], time) + 20], Y1[nearest_time(τ[1], time)]), (0., Y1[nearest_time(τ[2], time) - 1] - Y1[nearest_time(τ[1], time)]), color = :blue)
    text!(time[nearest_time(τ[1], time) + 20], (Y1[nearest_time(τ[1], time)] + Y1[nearest_time(τ[2], time) - 1]) / 2, text=L"x", color=:blue, fontsize=20.)

    # Grey dashed lines to illustrate the values of the process at maintenance dates
    lines!(time[[1, nearest_time(τ[2], time)]], [Y∞[nearest_time(τ[2], time) - 1] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(time[[nearest_time(τ[1], time), nearest_time(τ[2], time)]], [Y1[nearest_time(τ[1], time)] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[1]], linestyle=:dash, color=:magenta, linewidth=3.)
    maint_type2 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[2]], linestyle=:dash, color=:black, linewidth=3.)
    maint_type3 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[3]], linestyle=:dash, color=:green, linewidth=3.)

    axislegend(ax, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1 = 0.4", L"\rho_2 = 0.9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)
    ylims!(0., 25.)

    # Display the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\degradation_process_ARD1only.png", fig)


####MAINTENANCE EFFECT, FIRST EXAMPLE:ONE MAINTENANCE, ONLY ARD∞#####
fig = Figure()

    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "Time", ylabel = "Degradation Value")
    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    # Plot some data on each axis
    ARDinf_plot = lines!(ax, time, Y∞, color=:red, linewidth=1.)

    # ARD∞ arrows
    arrows!(ax, (time[nearest_time(τ[2], time) - 10], Y∞[nearest_time(τ[2], time) - 1]), (0., Y∞[nearest_time(τ[2], time)] - Y∞[nearest_time(τ[2], time) - 1]), color = :red)
    text!(time[nearest_time(τ[2], time) - 55], (Y∞[nearest_time(τ[2], time) - 1] + Y∞[nearest_time(τ[2], time)]) / 2 - 3., text=L"\rho_2 x", color=:red, fontsize=20.)
    arrows!(ax, (time[nearest_time(τ[1], time) + 10], 0.), (0., Y∞[nearest_time(τ[2], time) - 1]), color = :red)
    text!(time[nearest_time(τ[1], time) + 10], Y∞[nearest_time(τ[2], time) - 1] / 4, text=L"x", color=:red, fontsize=20.)

    # Grey dashed lines to illustrate the values of the process at maintenance dates
    lines!(time[[1, nearest_time(τ[2], time)]], [Y∞[nearest_time(τ[2], time) - 1] for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)
    lines!(time[[1, nearest_time(τ[2], time)]], [0. for _ in 1:2], linestyle=:dash, color=:grey, alpha=:.5)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[1]], linestyle=:dash, color=:magenta, linewidth=3.)
    maint_type2 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[2]], linestyle=:dash, color=:black, linewidth=3.)
    maint_type3 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[3]], linestyle=:dash, color=:green, linewidth=3.)

    axislegend(ax, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1 = 0.4", L"\rho_2 = 0.9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)
    ylims!(0., 25.)
    
    # Display the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\degradation_process_ARDinfonly.png", fig)



#ILLUSTRATES SAMPLING, ADD OF THE OBSERVATION POINTS, ONLY ARD1, PLUS SPECIAL OBSERVATION#####
fig = Figure()
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "", ylabel = "Degradation Value") 

    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    ARD1_degradation_process_plot = lines!(ax, time, Y1, linewidth=1., color=:blue)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[1]], linestyle=:dash, color=:magenta, linewidth=3.)
    maint_type2 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[2]], linestyle=:dash, color=:black, linewidth=3.)
    maint_type3 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[3]], linestyle=:dash, color=:green, linewidth=3.)

    I = [56, 102, 190, 270, 464, 498]
    observations = scatter!(ax, time[I], Y1[I], color=:black, markersize=13.)

    axislegend(ax, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1 = 0.4", L"\rho_2=0.9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)
    ylims!(0., 25.)
    
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\degradation_process_plus_observations_without_obs_notation.png", fig)



#ILLUSTRATES SAMPLING, ADD OF THE OBSERVATION POINTS, ONLY ARD1, PLUS SPECIAL OBSERVATION#####
fig = Figure()
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "", ylabel = "Degradation Value")

    axt = Axis(fig[1, 1], 
        xaxisposition=:bottom, xticks=([time[270]], [L"t_{\nu_3, 1}"]), xgridvisible=true, xticklabelsize=30.,
        yaxisposition=:right, yticks=([Y1[270]], [L"Y_{\nu_3, 1}"]), ygridvisible=true, yticklabelsize=30.)
    linkxaxes!(ax, axt)
    linkyaxes!(ax, axt)    

    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    ARD1_degradation_process_plot = lines!(ax, time, Y1, linewidth=1., color=:blue)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[1]], linestyle=:dash, color=:magenta, linewidth=3.)
    maint_type2 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[2]], linestyle=:dash, color=:black, linewidth=3.)
    maint_type3 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[3]], linestyle=:dash, color=:green, linewidth=3.)

    I = [56, 102, 190, 270, 464, 498]
    observations = scatter!(ax, time[I], Y1[I], color=:black, markersize=13.)

    axislegend(ax, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1 = 0.4", L"\rho_2=0.9", L"\rho_3 = 0.6"], position = :lt, labelsize=20.)
    ylims!(0., 25.)
    
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\degradation_process_plus_observations.png", fig)

#ILLUSTRATES SAMPLING, ADD OF THE OBSERVATION POINTS, ONLY ARD1, PLUS SPECIAL OBSERVATION#####
fig = Figure()
    ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false, xlabel = "Time", ylabel = "Degradation Value")

    axτ = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=20.,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax, axτ)
    linkyaxes!(ax, axτ)

    ARD1_degradation_process_plot = lines!(ax, time, Y1, linewidth=1., color=:blue)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[1]], linestyle=:dash, color=:magenta, linewidth=3.)
    maint_type2 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[2]], linestyle=:dash, color=:black, linewidth=3.)
    maint_type3 = vlines!(ax, [τ[i] for i in eachindex(τ) if maintenances[i, "TYPE"] == maint_types[3]], linestyle=:dash, color=:green, linewidth=3.)

    I = [56, 102, 190, 270, 464, 498]
    observations = scatter!(ax, time[I], Y1[I], color=:black, markersize=13.)

    lines!(time[[102, 190]], Y1[[190, 190]], linestyle=:dot, color=:grey, alpha=:.5)
    bracket!(ax, (time[102], Y1[102]), (time[102], Y1[190]), color = :red, rotation = 0)
    text!(time[30], (Y∞[190] + Y1[102])/2 - 0.2, text=L"\Delta Y_{\nu_2, 1}", color=:red, fontsize=18.)

    lines!(time[[270, 190]], Y1[[190, 190]], linestyle=:dot, color=:grey, alpha=:.5)
    bracket!(ax, (time[270], Y1[190]), (time[270], Y1[270]), color = RGBAf(52/255, 119/255, 37/255, 1.), rotation = 0)
    text!(time[300], (Y∞[190] + Y1[270])/2 - 0.6, text=L"Z_{\nu_2}", color = RGBAf(52/255, 119/255, 37/255, 1.), fontsize=18.)

    #    axislegend(ax, merge = true, [ARD1_degradation_process_plot, maint_type1, maint_type2, maint_type3, observations], [L"ARD_1", L"\rho_1 = 0.4", L"\rho_2=1.", L"\rho_3 = 0.6", L"Observations"], position = :lt, labelsize=20.)
    ylims!(0., 25.)
    
    display(fig)
save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\degradation_process_plus_observations_illustration_incrment_plus_jump.png", fig)



#####gaps between maintenances#####
#ARD1 parameters
μ = 2.
σ = sqrt(5)
ρ = [.5]

#create the wienerARD1 instance
wienerARD1 = WienerARD1(["P"], μ, σ, ρ)

#degradationdata parameters
k = 3#nb of maintenances
nj = 3#nb of obseravtions between each maintenances
T = (k+1)*nj
τ = zeros(k)
for i in 1:k
    τ[i] = i*(nj+1)
end
τ
τ_types = ["P", "P", "P"]
T = τ[end] + nj
Δt = 1.

#create the degradationdata instance
degradationdata = CreateData2(τ, Δt, T, τ_types)
Random.seed!(11)
rand!(wienerARD1, degradationdata, 0)

#let us cancel some data to make it more asymetric
maint = degradationdata.maintenances
deg = degradationdata.degradations
deg = deg[Not(5, 11, 12, 13, 17), :]
degradationdata = DegradationData(maint, deg)
deg_obs = filter(row -> row.TYPE in ["Between"], degradationdata.degradations)
degradationdata_obs = DegradationData(maint, deg_obs)

#creates two separated sets of dataframes (so it is easy to show)
deg_unobs2 = filter(row -> row.NB_MAINTENANCES >= 1 && (row.NB_MAINTENANCES != 2 || row.TYPE == "Before" || row.TYPE == "After"), deg)
deg_obs2 = filter(row -> row.NB_MAINTENANCES >=1, deg_obs)



#####PLOT TIME#####
#first plot for observed increments
fig = Figure()

    #axes
    ax2 = Axis(fig[1, 1], xlabel = "Observation Date", xticks = (vcat([0.], deg_obs2.DATE), vcat([L"t_{0, 0} = 0 ⋯"], [L"t_{\nu_i,%$i}" for i in 1:length(deg_obs2.DATE) - 2], [L"t_{\nu_{i+1},%$i}" for i in 1:2])), xgridvisible=false, xticklabelsize=20.,  
        yticks = (deg_obs2.VALUE, ["hello" for i in eachindex(deg_obs2.VALUE)]), yticklabelsvisible = false, ygridvisible=true)
    ax3 = Axis(fig[1, 1], xticklabelcolor = :black, xaxisposition = :top, xlabel = "Maintenance Date", xticks = (τ, [L"\tau_{\nu_{i}-1}", L"\tau_{\nu_i}", L"\tau_{\nu_i+1} = \tau_{\nu_{i+1}-1}"]), xticklabelpad = 10.0, yticks = ([], []), yticklabelsvisible = false, xgridvisible = false, xticklabelsize=20.)
    linkxaxes!(ax2, ax3)
    linkyaxes!(ax2, ax3)

    #to show maintenances
    for t in τ
        vlines!(ax2, t, color=:black, linestyle=:dash)
    end

    #trajectory
    lines!(ax2, deg_obs2.DATE, deg_obs2.VALUE, linewidth = 2, color=:black)
    l1 = bracket!(ax2, (deg_obs2.DATE[1], deg_obs2.VALUE[1]), (deg_obs2.DATE[1], deg_obs2.VALUE[2]), color=:red, text=L"\Delta y_{\nu_i, 2}")
    l2 = bracket!(ax2, (deg_obs2.DATE[2], deg_obs2.VALUE[2]), (deg_obs2.DATE[2], deg_obs2.VALUE[3]), color=:red, text=L"\Delta y_{\nu_i, 3}")
    l3 = bracket!(ax2, (deg_obs2.DATE[3], deg_obs2.VALUE[3]), (deg_obs2.DATE[3], deg_obs2.VALUE[4]), color=:green, text=L"z_{\nu_i}")
    l4 = bracket!(ax2, (deg_obs2.DATE[4], deg_obs2.VALUE[4]), (deg_obs2.DATE[4], deg_obs2.VALUE[5]), color=:red, text=L"\Delta y_{\nu_{i+1}, 2}")
    scatter!(ax2, deg_obs2.DATE, deg_obs2.VALUE, color=:black)

    display(fig)

    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\observed_increments.png", fig)

#second plot for unobserved increments
fig = Figure()

    #axes
    ax2 = Axis(fig[1, 1], xlabel = "Observation Date", xticks = (vcat([0.], deg_obs2.DATE), vcat([L"t_{0, 0} = 0 ⋯"], [L"t_{\nu_i,%$i}" for i in 1:length(deg_obs2.DATE) - 2], [L"t_{\nu_{i+1},%$i}" for i in 1:2])), xgridvisible=false, xticklabelsize=20.,
        yticks = (deg_unobs2.VALUE[3:8], ["hello" for i in eachindex(deg_unobs2.VALUE[3:8])]), yticklabelsvisible = false, ygridvisible=true)
    ax3 = Axis(fig[1, 1], xticklabelcolor = :black, xaxisposition = :top, xlabel = "Maintenance Date", xticks = (τ, [L"\tau_{\nu_{i}-1}", L"\tau_{\nu_i}", L"\tau_{\nu_i+1} = \tau_{\nu_{i+1}-1}"]), xticklabelpad = 10.0, yticks = ([], []), yticklabelsvisible = false, xgridvisible = false, xticklabelsize=20.)
    linkxaxes!(ax2, ax3)
    linkyaxes!(ax2, ax3)

    #to show maintenances
    for t in τ
        vlines!(ax2, t, color=:black, linestyle=:dash)
    end

    #trajectory
    #l0 = lines!(ax2, deg_obs2.DATE, deg_obs2.VALUE, color = :black, linewidth = 2.)
    l1 = lines!(ax2, deg_unobs2.DATE, deg_unobs2.VALUE, 
        color = :blue, linestyle=:dash, linewidth = 1.)
    l2 = bracket!(ax2, (deg_unobs2.DATE[3], deg_unobs2.VALUE[3]), (deg_unobs2.DATE[3], deg_unobs2.VALUE[4]), 
        color = :orange, text=L"\Delta Y_{\nu_i, n_{\nu_i}+1}")
    l3 = bracket!(ax2, (deg_unobs2.DATE[5], deg_unobs2.VALUE[5]), (deg_unobs2.DATE[4], deg_unobs2.VALUE[4]), 
        color = :lightgreen, text=L"\mathcal{Z}_{\nu_i}^\ast")
    l4 = bracket!(ax2, ((deg_unobs2.DATE[5] + deg_unobs2.DATE[6]) / 2, deg_unobs2.VALUE[5]), ((deg_unobs2.DATE[5] + deg_unobs2.DATE[6]) / 2, deg_unobs2.VALUE[6]), 
        color = :orange, text=L"\Delta Y_{\nu_{i+1}-1, 1}")
    l5 = bracket!(ax2, (deg_unobs2.DATE[7], deg_unobs2.VALUE[7]), (deg_unobs2.DATE[6], deg_unobs2.VALUE[6]), 
        color = :lightgreen, text=L"\mathcal{Z}_{\nu_{i}+1}^\ast")
    l6 = bracket!(ax2, (deg_unobs2.DATE[8], deg_unobs2.VALUE[8]), (deg_unobs2.DATE[8], deg_unobs2.VALUE[7]), 
        color = :orange, text=L"\Delta Y_{\nu_{i+1}, 1}")
    scatter!(ax2, deg_unobs2.DATE, deg_unobs2.VALUE, color=:blue, markersize=6.)
    scatter!(ax2, deg_obs2.DATE, deg_obs2.VALUE, color=:black)

    display(fig)

    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\unobserved_increments.png", fig)



#degradationdata parameters
k = 3#nb of maintenances
nj = 3#nb of obseravtions between each maintenances
T = (k+1)*nj
τ = zeros(k)
for i in 1:k
    τ[i] = i*(nj+1)
end
τ
τ_types = ["P", "P", "P"]
T = τ[end] + nj
Δt = 1.

#create the degradationdata instance
degradationdata = CreateData2(τ, Δt, T, τ_types)
Random.seed!(12)
rand!(wienerARD1, degradationdata, 0)

#let us cancel some data to make it more asymetric
maint = degradationdata.maintenances
deg = degradationdata.degradations
deg = deg[Not(5, 9, 11, 12, 17), :]
degradationdata = DegradationData(maint, deg)
deg_obs = filter(row -> row.TYPE in ["Between"], degradationdata.degradations)
degradationdata_obs = DegradationData(maint, deg_obs)

#creates two separated sets of dataframes (so it is easy to show)
deg_unobs2 = filter(row -> row.NB_MAINTENANCES >= 1, deg)
deg_obs2 = filter(row -> row.NB_MAINTENANCES >=1, deg_obs)
degradationdata_unobs2 = DegradationData(maint, deg_unobs2)



#####PLOT TIME#####
#first plot for observed increments
fig = Figure()

    #axes
    ax2 = Axis(fig[1, 1], xlabel = "Observation Date", xticks = (vcat([0.], deg_obs2.DATE), vcat([L"t_{0, 0} = 0 ⋯"], [L"t_{\nu_i,%$i}" for i in 1:length(deg_obs2.DATE) - 2], [L"t_{\nu_{i+1},%$i}" for i in 1:2])), xgridvisible=false, xticklabelsize=20.,  
        yticks = (deg_obs2.VALUE, ["hello" for i in eachindex(deg_obs2.VALUE)]), yticklabelsvisible = false, ygridvisible=true)
    ax3 = Axis(fig[1, 1], xticklabelcolor = :black, xaxisposition = :top, xlabel = "Maintenance Date", xticks = (τ, [L"\tau_{\nu_{i}-1}", L"\tau_{\nu_i}", L"\tau_{\nu_i+1} = \tau_{\nu_{i+1}-1}"]), xticklabelpad = 10.0, yticks = ([], []), yticklabelsvisible = false, xgridvisible = false, xticklabelsize=20.)
    linkxaxes!(ax2, ax3)
    linkyaxes!(ax2, ax3)

    #to show maintenances
    for t in τ
        vlines!(ax2, t, color=:black, linestyle=:dash)
    end

    #trajectory
    lines!(ax2, deg_unobs2.DATE, deg_unobs2.VALUE, linewidth = 2, color=:black, linestyle=:dash)
    #scatter!(ax2, deg_unobs2.DATE, deg_unobs2.VALUE, color=:black)
    plot!(ax2, degradationdata_unobs2, linecolor=:black)

    display(fig)

    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\non_invertible_case.png", fig)


#third plot to show what happens if ν_1 = 1 or not
#degradationdata parameters
k = 5#nb of maintenances
nj = 3#nb of obseravtions between each maintenances
T = (k+1)*nj
τ = zeros(k)
for i in 1:k
    τ[i] = i*(nj+1)
end
τ
τ_types = ["P", "P", "P", "P", "P"]
T = τ[end] + nj
Δt = 1.

#create the degradationdata instance
degradationdata = CreateData2(τ, Δt, T, τ_types)
Random.seed!(4)
rand!(wienerARD1, degradationdata, 0)

#We want one more "Before" observation so we need to create 8 maintenances and cancel all data linked to the eight-th maintenance
maint = degradationdata.maintenances[1:2, :]
deg = filter(row -> row.NB_MAINTENANCES <= 2, degradationdata.degradations)

#let us cancel some data to make it more asymetric
deg = deg[Not([2, 3, 4, 7, 13, 14]), :]
degradationdata = DegradationData(maint, deg)
deg_obs = filter(row -> row.TYPE in ["Between"] || (row.TYPE == "After" && row.NB_MAINTENANCES in [5]) || (row.TYPE == "Before" && row.NB_MAINTENANCES in [0]), degradationdata.degradations)
degradationdata_obs = DegradationData(maint, deg_obs)

#creates two separated sets of dataframes (so it is easy to show)
deg_without0 = filter(row -> row.NB_MAINTENANCES != 0, deg_obs)

custom_theme = Theme(Axis = (ygridvisible = false,),)
set_theme!(custom_theme)

fig = Figure()

    #axes
    ax1 = Axis(fig[1, 1], xticks = (deg_obs.DATE, [L"t_{%$i}" for i in 1:length(deg_obs.DATE)]), xgridvisible = false, xticklabelsize=20., ylabel = "Degradation Value")
    ax2 = Axis(fig[1, 2], xticks = (deg_without0.DATE, [L"t_{%$i}" for i in 1:length(deg_without0.DATE)]), xgridvisible = false, xticklabelsize=20., ylabel = "Degradation Value")
    ax3 = Axis(fig[1, 1], xticklabelcolor = :black, xaxisposition = :top, xlabel = L"\nu_1 = 1", xticks = (τ[1:nrow(maint)], [L"\tau_{%$i}" for i in 1:nrow(maint)]), xticklabelsize=20., yticks = ([], []), yticklabelsvisible = false, ygridvisible = false)
    ax4 = Axis(fig[1, 2], xticklabelcolor = :black, xaxisposition = :top, xlabel = L"\nu_1 \in \mathbb{N} \setminus \{ 0, 1 \}", xticks = (τ[1:nrow(maint)], [L"\tau_{%$i}" for i in 1:nrow(maint)]), xticklabelsize=20., yticks = ([], []), yticklabelsvisible = false, ygridvisible = false)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    linkxaxes!(ax2, ax4)
    linkxaxes!(ax1, ax3)

    #to show maintenances
    for t in τ[1:nrow(maint)]
        vlines!(ax1, t, color=:black, linestyle=:dash)
        vlines!(ax2, t, color=:black, linestyle=:dash)
    end

    #trajectory
    lines!(ax1, deg_obs.DATE, deg_obs.VALUE, color = :black, linewidth = 2)
    lines!(ax2, deg_without0.DATE, deg_without0.VALUE, color = :black, linewidth = 2)
    scatter!(ax1, deg_obs.DATE, deg_obs.VALUE, color=:black)
    scatter!(ax2, deg_without0.DATE, deg_without0.VALUE, color=:black)

    #premier truc
    lines!(ax1, [0., deg_obs.DATE[1]], [0., deg_obs.VALUE[1]], linestyle=:dash, color = :red, linewidth = 2)
    bracket!(ax1, (deg_obs.DATE[1], deg_obs.VALUE[1]), (deg_obs.DATE[1], 0.), color = :red, text=L"\Delta Y_{\nu_1, 1} = \Delta Y_{1, 1}")
    lines!(ax2, [0., deg_without0.DATE[1]], [0., deg_without0.VALUE[1]], linestyle=:dash, color = :green, linewidth = 2)
    bracket!(ax2, (deg_without0.DATE[1], deg_without0.VALUE[1]), (deg_without0.DATE[1], 0.), color = :green, text=L"Z_{\nu_0} = Z_0")

    display(fig)

    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\nu_1egal1.png", fig)