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
mward1 = ARD.MWARD1([2., 2.], [.4 .2; .2 .4], Dict("P" => [0.5, 0.2], "C" => [0.8, 0.5], "T" => [.9, .9]))
mward∞ = ARD.MWARD∞([2., 2.], [.4 .2; .2 .4], Dict("P" => [0.5, 0.2], "C" => [0.8, 0.5], "T" => [.9, .9]))

μ = [2., 2.]
Σ = [.4 .2; .2 .4]
ρ = Dict((:ind1, "P") => ARD.Efficiency(.5, ARD.ARD1()), (:ind2, "P") => ARD.Efficiency(.2, ARD.ARDinf()), (:ind1, "C") => ARD.Efficiency(.8, ARD.ARD1()), (:ind2, "C") => ARD.Efficiency(.5, ARD.ARDinf()), (:ind1, "T") => ARD.Efficiency(.9, ARD.ARD1()), (:ind2, "T") => ARD.Efficiency(.9, ARD.ARDinf()))
mvw = ARD.MvWienerAR(μ, Σ, ρ)

# Define the inspection dates
T = 15.
inspection_dates = convert(Vector{Float64}, range(0, T, 500))

# Maintenance dates
τ = T .* [.17, .4, .5, .8]
maintenances = DataFrame(DATE=τ, TYPE=["P", "C", "P", "T"])

# Generation of the random processes trajectories
nb_seed = 3
Random.seed!(nb_seed)
X = ARD.mw_rand(mvw, inspection_dates)
Random.seed!(nb_seed)
Y = ARD.rand(mvw, inspection_dates, maintenances)
Y1 = Y[1, :]
Y∞ = Y[2, :]

# Global variable for label size
lbsize = 40.
lbsize2 = 30.

##### BIVARIATE UNMAINTAINED PROCESS #####
# Two figures, one for each characteristic
fig = Figure(resolution = (1920, 1080))

    # Axis definition
    ax1 = Axis(fig[1, 1], ylabel = "Characteristic (1)", ylabelsize=lbsize2, xgridvisible=false)
    ax2 = Axis(fig[2, 1], ylabel = "Characteristic (2)", ylabelsize=lbsize2, xlabel = "Time", xlabelsize=lbsize2, xgridvisible=false)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, X[1, :], color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, X[2, :], color=:red, linewidth=1.)

    # Drift directing curve
    drift_char1 = lines!(ax1, inspection_dates, mvw.drift[1] .* inspection_dates, color=:blue, alpha=0.5)
    drift_char2 = lines!(ax2, inspection_dates, mvw.drift[2] .* inspection_dates, color=:red, alpha=0.5)

    ylims!(ax1, -2.5, 37.)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_unmaintained_2fig.png", fig)



# One figure with both characteristics
fig = Figure(resolution = (1920, 1080))

    # Axis definition
    ax1 = Axis(fig[1, 1], ylabel = "Degradation level", ylabelsize=lbsize2, xlabel = "Time", xlabelsize=lbsize2, xgridvisible=false)

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, X[1, :], color=:blue, linewidth=1.)
    char2 = lines!(ax1, inspection_dates, X[2, :], color=:red, linewidth=1.)

    # Drift directing curve
    drift_char1 = lines!(ax1, inspection_dates, mvw.drift[1] .* inspection_dates, color=:blue, alpha=0.5)
    drift_char2 = lines!(ax1, inspection_dates, mvw.drift[2] .* inspection_dates, color=:red, alpha=0.5)

    ylims!(ax1, -2.5, 40.)

    # Display and save the figure
    display(fig)
save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_unmaintained_1fig.png", fig)



μ_list = [[1., 2.], [2., 2.], [2., 2.], [2., 2.]]
Σ_list = [[.4 .2; .2 .4], [.4 .37; .37 .4], [2. .2; .2 .3], [.4 0.; 0. .4]]

fig = Figure(resolution = (1920, 1080))
vect_ax = [Axis(fig[1, 1], xgridvisible=false, xlabel=L"(a)", xlabelsize=lbsize2), Axis(fig[1, 2], xgridvisible=false, xlabel=L"(b)", xlabelsize=lbsize2), Axis(fig[2, 1], xgridvisible=false, xlabel=L"(c)", xlabelsize=lbsize2), Axis(fig[2, 2], xgridvisible=false, xlabel=L"(d)", xlabelsize=lbsize2)]
for i in eachindex(μ_list)
    mvw = ARD.MvWienerAR(μ_list[i], Σ_list[i], ρ)
    Random.seed!(nb_seed)
    X = ARD.mw_rand(mvw, inspection_dates)
    char1 = lines!(vect_ax[i], inspection_dates, X[1, :], color=:blue, linewidth=1.)
    char2 = lines!(vect_ax[i], inspection_dates, X[2, :], color=:red, linewidth=1.)
    drift_char1 = lines!(vect_ax[i], inspection_dates, mvw.drift[1] .* inspection_dates, color=:blue, alpha=0.5)
    drift_char2 = lines!(vect_ax[i], inspection_dates, mvw.drift[2] .* inspection_dates, color=:red, alpha=0.5)
    ylims!(vect_ax[i], -2.5, 40.)
end
display(fig)

save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_unmaintained_4mw.png", fig)



##### BIVARIATE ARD1 PROCESS WITH FIVE MAINTENANCE ACTIONS, MAIN FIGURE ######
# Inspection dates indices
I = [18, 56, 102, 270, 410, 464, 500]
I_char1 = I[[3, 4, 5, 7]]
I_char2 = I[[1, 2, 6, 7]]

# ARD1 figure
fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1], xticks=(inspection_dates[I_char1], [L"t_{\nu_1^{(1)}, \phi_{1, 1}^{(1)}}", L"t_{\nu_2^{(1)}, \phi_{2, 1}^{(1)}}", L"t_{\nu_3^{(1)}, \phi_{3, 1}^{(1)}}", L"t_{\nu_3^{(1)}, \phi_{3, 2}^{(1)}}"]), ylabel = L"\textrm{Characteristic } (1): ARD_1", ylabelsize=lbsize2, xgridvisible=false, xticklabelsize=lbsize, ygridvisible=false)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Inspection Date", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1], xticks=(inspection_dates[I_char2], [L"t_{\nu_1^{(2)}, \phi_{1, 1}^{(2)}}", L"t_{\nu_1^{(2)}, \phi_{1, 2}^{(2)}}", L"t_{\nu_2^{(2)}, \phi_{2, 1}^{(2)}}", L"t_{\nu_2^{(2)}, \phi_{2, 2}^{(2)}}"]), ylabel = L"\textrm{Characteristic } (2): ARD_\infty", ylabelsize=lbsize2, xgridvisible=false, xticklabelsize=lbsize, ygridvisible=false)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, Y∞, color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax2, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y∞[I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    vlines!(ax2, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=lbsize2)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=lbsize2)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure.png", fig)


# No ticks for the axis of fig 1 and 3
# ARD1 figure
fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (1): ARD_1", ygridvisible=false, ylabelsize=lbsize2)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Inspection Date", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (2): ARD_\infty", ygridvisible=false, ylabelsize=lbsize2)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, Y∞, color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax2, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y∞[I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    vlines!(ax2, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=lbsize2)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=lbsize2)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure_without_phi.png", fig)


# illustration of the maintenance effects with arrows
fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (1): ARD_1", ygridvisible=false, ylabelsize=lbsize2)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Inspection Date", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (2): ARD_\infty", ygridvisible=false, ylabelsize=lbsize2)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, Y∞, color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax2, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y∞[I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    vlines!(ax2, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=lbsize2)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=lbsize2)

    # ARD1 arrows
    arrows2d!(ax1, (inspection_dates[nearest_time(τ[4], inspection_dates)], Y1[nearest_time(τ[4], inspection_dates) - 1]), (0., Y1[nearest_time(τ[4], inspection_dates)] - Y1[nearest_time(τ[4], inspection_dates) - 1]), color = :blue, tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8)
    text!(ax1, inspection_dates[nearest_time(τ[4], inspection_dates) - 70], (Y1[nearest_time(τ[4], inspection_dates) - 1] + Y1[nearest_time(τ[4], inspection_dates)]) / 2 - 1.5, text=L"\mathcal{Z}_4^{(1)} = - \rho_2 x", color=:blue, fontsize=lbsize)
    arrows2d!(ax1, (inspection_dates[nearest_time(τ[3], inspection_dates)], Y1[nearest_time(τ[3], inspection_dates)]), (0., Y1[nearest_time(τ[4], inspection_dates) - 1] - Y1[nearest_time(τ[3], inspection_dates)]), color = :blue)
    text!(ax1, inspection_dates[nearest_time(τ[3], inspection_dates) + 10], (Y1[nearest_time(τ[3], inspection_dates)] + Y1[nearest_time(τ[4], inspection_dates) - 1]) / 2, text=L"x", color=:blue, fontsize=lbsize)
    hlines!(ax1, Y1[nearest_time(τ[4], inspection_dates) - 1], color=:blue, linestyle=:dash, linewidth=1.)
    hlines!(ax1, Y1[nearest_time(τ[3], inspection_dates)], color=:blue, linestyle=:dash, linewidth=1.)

    # ARD∞ arrows
    arrows2d!(ax2, (inspection_dates[nearest_time(τ[4], inspection_dates)], Y∞[nearest_time(τ[4], inspection_dates) - 1]), (0., Y∞[nearest_time(τ[4], inspection_dates)] - Y∞[nearest_time(τ[4], inspection_dates) - 1]), color = :red, tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8)
    text!(ax2, inspection_dates[nearest_time(τ[4], inspection_dates) - 70], (Y∞[nearest_time(τ[4], inspection_dates) - 1] + Y∞[nearest_time(τ[4], inspection_dates)]) / 2 - 1, text=L"\mathcal{Z}_4^{(2)} = - \rho_2 x", color=:red, fontsize=lbsize)
    arrows2d!(ax2, (inspection_dates[nearest_time(τ[3], inspection_dates)], 0.), (0., Y∞[nearest_time(τ[4], inspection_dates) - 1]), color = :red)
    text!(ax2, inspection_dates[nearest_time(τ[3], inspection_dates) + 10], (Y∞[nearest_time(τ[3], inspection_dates) - 1] + Y∞[nearest_time(τ[4], inspection_dates) - 1]) / 2, text=L"x", color=:red, fontsize=lbsize)
    hlines!(ax2, Y∞[nearest_time(τ[4], inspection_dates) - 1], color=:red, linestyle=:dash, linewidth=1.)
    hlines!(ax2, 0., color=:red, linestyle=:dash, linewidth=1.)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure_ARD_illustration.png", fig)


# Text for number of inspections, without the φ in the x-axis labels
# ARD1 figure
fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (1): ARD_1", ygridvisible=false, ylabelsize=lbsize2)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel=L"\textrm{Maintenance Date}-\mathcal{K} = 5", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=lbsize, xlabelsize=lbsize,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Inspection Date", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (2): ARD_\infty", ygridvisible=false, ylabelsize=lbsize2)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, Y∞, color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax2, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y∞[I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    vlines!(ax2, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=lbsize2)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=lbsize2)

    # Add text for number of observations
    text!(ax1, (0. + τ[1]) / 2 - 1., 15., text=L"\mathcal{N}_1 = 2", color=:black, fontsize=lbsize)
    text!(ax1, (τ[1] + τ[2]) / 2 - .5, 15., text=L"\mathcal{N}_2 = 1", color=:black, fontsize=lbsize)
    text!(ax1, (τ[2] + τ[3]) / 2 - .5, 15., text=L"\mathcal{N}_3 = 0", color=:black, fontsize=lbsize)
    text!(ax1, (τ[3] + τ[4]) / 2 - .5, 15., text=L"\mathcal{N}_4 = 1", color=:black, fontsize=lbsize)
    text!(ax1, (τ[4] + T) / 2 - 1., 15., text=L"\mathcal{N}_5 = 3", color=:black, fontsize=lbsize)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure_virtual_numerotation.png", fig)


# Text for number of inspections, without the φ in the x-axis labels
# ARD1 figure
fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (1): ARD_1", ygridvisible=false, ylabelsize=lbsize2)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel=L"\textrm{Maintenance Date}-\mathcal{K} = 5", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=lbsize, xlabelsize=lbsize,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Inspection Date", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1],
        xticks=([], []), xgridvisible=false,
        ylabel = L"\textrm{Characteristic } (2): ARD_\infty", ygridvisible=false, ylabelsize=lbsize2)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, Y∞, color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax2, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y∞[I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    vlines!(ax2, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=lbsize2)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=lbsize2)

    # Add text for number of observations for the first characteristic
    text!(ax1, [((τ[1] + τ[2]) / 2 - .5, 15.), ((τ[1] + τ[2]) / 2 - .5, 25.)], text=[L"N_1^{(1)} = 1", L"\nu_1^{(1)} = 2"], color=:black, fontsize=lbsize)
    text!(ax1, [((τ[3] + τ[4]) / 2 - .5, 15.), ((τ[3] + τ[4]) / 2 - .5, 25.)], text=[L"N_2^{(1)} = 1", L"\nu_2^{(1)} = 4"], color=:black, fontsize=lbsize)
    text!(ax1, [((τ[4] + T) / 2 - .5, 15.), ((τ[4] + T) / 2 - .5, 25.)], text=[L"N_3^{(1)} = 2", L"\nu_3^{(1)} = 5"], color=:black, fontsize=lbsize)

    # Add text for number of observations for the second characteristic
    text!(ax2, [((0. + τ[1]) / 2 - .5, 15.), ((0. + τ[1]) / 2 - .5, 25.)], text=[L"N_1^{(2)} = 2", L"\nu_1^{(2)} = 1"], color=:black, fontsize=lbsize)
    text!(ax2, [((τ[4] + T) / 2 - .5, 15.), ((τ[4] + T) / 2 - .5, 25.)], text=[L"N_2^{(2)} = 2", L"\nu_2^{(2)} = 5"], color=:black, fontsize=lbsize)

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure_observed_numerotation.png", fig)



# One figure with both characteristics
fig = Figure(resolution = (1920, 1080))

    # Axis definition
    ax1 = Axis(fig[1, 1], ylabel = "Degradation level", ylabelsize=lbsize2, xlabel = "Time", xlabelsize=lbsize2, xgridvisible=false)

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, X[1, :], color=:blue, linewidth=1.)
    char2 = lines!(ax1, inspection_dates, X[2, :], color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 40.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], X[1, I_char1], color=:black, markersize=20.)
    observations_char1 = scatter!(ax1, inspection_dates[I_char2], X[2, I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Add text for number of observations
    text!(ax1, (0. + τ[1]) / 2 - 1., 15., text=L"\mathcal{N}_1 = 2", color=:black, fontsize=lbsize)
    text!(ax1, (τ[1] + τ[2]) / 2 - .5, 15., text=L"\mathcal{N}_2 = 1", color=:black, fontsize=lbsize)
    text!(ax1, (τ[2] + τ[3]) / 2 - .5, 15., text=L"\mathcal{N}_3 = 0", color=:black, fontsize=lbsize)
    text!(ax1, (τ[3] + τ[4]) / 2 - .5, 15., text=L"\mathcal{N}_4 = 1", color=:black, fontsize=lbsize)
    text!(ax1, (τ[4] + T) / 2 - 1., 15., text=L"\mathcal{N}_5 = 3", color=:black, fontsize=lbsize)

    #all virtual increments
    ps=Point2f.([(t, 0) for t in sort(vcat(0., inspection_dates[I], τ))])
    text_virtual_var = [L"\Delta X_{1, 1}", L"\Delta X_{1, 2}", L"\Delta X_{1, 3}", L"\Delta X_{2, 1}", L"\Delta X_{2, 2}", L"\Delta X_{3, 1}", L"\Delta X_{4, 1}", L"\Delta X_{4, 2}", L"\Delta X_{5, 1}", L"\Delta X_{5, 2}", L"\Delta X_{5, 3}"]
    for i in 1:length(ps)-1
        arrows2d!(ax1, ps[i], ps[i+1],
        argmode = :endpoint,
        tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8
        )
        text!(ax1, (ps[i][1] + ps[i+1][1]) / 2 - .5, -2, text=text_virtual_var[i], color=:black, fontsize=lbsize - 5)
    end

    # virtual increments for the sole last observations of the second characteristic
    #ps=Point2f.([(t, 0) for t in inspection_dates[I][end-2:end]])
    #for i in 1:length(ps)-1
    #    arrows2d!(ax1, ps[i], ps[i+1],
    #    argmode = :endpoint,
    #    tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8
    #    )
    #    text!(ax1, (ps[i][1] + ps[i+1][1]) / 2 - .5, -2, text=text_virtual_var[i], color=:black, fontsize=lbsize - 5)
    #end
I
    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_unmaintained_virtual_observations.png", fig)




# One figure with both characteristics
fig = Figure(resolution = (1920, 1080))

    # Axis definition
    ax1 = Axis(fig[1, 1], ylabel = "Degradation level", ylabelsize=lbsize2, xlabel = "Time", xlabelsize=lbsize2, xgridvisible=false)

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, X[1, :], color=:blue, linewidth=1., alpha=.6)
    char1ARD = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 40.)

    # Add observation points with red color and bigger markersize
    all_observations_char1 = scatter!(ax1, inspection_dates[I], X[1, I], color=:grey, markersize=20.)
    observations_char1_ARD = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Add text for number of observations
    text!(ax1, (0. + τ[1]) / 2 - 1., 15., text=L"\mathcal{N}_1 = 2", color=:black, fontsize=lbsize)
    text!(ax1, (τ[1] + τ[2]) / 2 - .5, 15., text=L"\mathcal{N}_2 = 1", color=:black, fontsize=lbsize)
    text!(ax1, (τ[2] + τ[3]) / 2 - .5, 15., text=L"\mathcal{N}_3 = 0", color=:black, fontsize=lbsize)
    text!(ax1, (τ[3] + τ[4]) / 2 - .5, 15., text=L"\mathcal{N}_4 = 1", color=:black, fontsize=lbsize)
    text!(ax1, (τ[4] + T) / 2 - 1., 15., text=L"\mathcal{N}_5 = 3", color=:black, fontsize=lbsize)

    #all virtual increments
    ps=Point2f.([(t, 30.) for t in sort(vcat(0., inspection_dates[I], τ))])
    text_virtual_var = [L"\Delta X_{1, 1}", L"\Delta X_{1, 2}", L"\Delta X_{1, 3}", L"\Delta X_{2, 1}", L"\Delta X_{2, 2}", L"\Delta X_{3, 1}", L"\Delta X_{4, 1}", L"\Delta X_{4, 2}", L"\Delta X_{5, 1}", L"\Delta X_{5, 2}", L"\Delta X_{5, 3}"]
    for i in 1:length(ps)-1
        arrows2d!(ax1, ps[i], ps[i+1],
        argmode = :endpoint,
        tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8
        )
        text!(ax1, (ps[i][1] + ps[i+1][1]) / 2 - .5, 30., text=text_virtual_var[i], color=:black, fontsize=lbsize - 5)
    end

    # observed increments for the sole last observations of the first characteristic
    ps1=Point2f.([(t, 0.) for t in vcat(0., inspection_dates[I_char1])])
    text_obs_var = [L"Z_0^{(1)}", L"Z_1^{(1)}", L"Z_2^{(1)}", L"\Delta Y_{3, 2}"]
    for i in 1:length(ps1)-1
        arrows2d!(ax1, ps1[i], ps1[i+1],
        argmode = :endpoint,
        tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8
        )
        text!(ax1, (ps1[i][1] + ps1[i+1][1]) / 2 - .5, 0., text=text_obs_var[i], color=:black, fontsize=lbsize - 5)
    end
I
    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_unmaintained_observed_increments_matrix.png", fig)




fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1], xticks=(inspection_dates[I_char1], [L"t_{\nu_1^{(1)}, \phi_{1, 1}^{(1)}}", L"t_{\nu_2^{(1)}, \phi_{2, 1}^{(1)}}", L"t_{\nu_3^{(1)}, \phi_{3, 1}^{(1)}}", L"t_{\nu_3^{(1)}, \phi_{3, 2}^{(1)}}"]), ylabel = L"\textrm{Characteristic } (1): ARD_1", ylabelsize=lbsize2, xgridvisible=false, xticklabelsize=lbsize, ygridvisible=false)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Inspection Date", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1], xticks=(inspection_dates[I_char2], [L"t_{\nu_1^{(2)}, \phi_{1, 1}^{(2)}}", L"t_{\nu_1^{(2)}, \phi_{1, 2}^{(2)}}", L"t_{\nu_2^{(2)}, \phi_{2, 1}^{(2)}}", L"t_{\nu_2^{(2)}, \phi_{2, 2}^{(2)}}"]), ylabel = L"\textrm{Characteristic } (2): ARD_\infty", ylabelsize=lbsize2, xgridvisible=false, xticklabelsize=lbsize, ygridvisible=false)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, Y∞, color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax2, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y∞[I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    vlines!(ax2, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=lbsize2)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=lbsize2)

    # virtual increments for the sole last observations of the second characteristic
    ps=Point2f.([(t, 0) for t in inspection_dates[I][end-2:end]])
    for i in 1:length(ps)-1
        arrows2d!(ax1, ps[i], ps[i+1],
        argmode = :endpoint,
        tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8
        )
        text!(ax1, (ps[i][1] + ps[i+1][1]) / 2 - .5, 0., text=text_virtual_var[i], color=:black, fontsize=lbsize - 5)
    end

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure_observed_increments_decompostition.png", fig)


fig = Figure(resolution = (1920, 1080))

    # Axis definition for the first figure
    ax1 = Axis(fig[1, 1], xticks=(inspection_dates[I_char1], [L"t_{\nu_1^{(1)}, \phi_{1, 1}^{(1)}}", L"t_{\nu_2^{(1)}, \phi_{2, 1}^{(1)}}", L"t_{\nu_3^{(1)}, \phi_{3, 1}^{(1)}}", L"t_{\nu_3^{(1)}, \phi_{3, 2}^{(1)}}"]), ylabel = L"\textrm{Characteristic } (1): ARD_1", ylabelsize=lbsize2, xgridvisible=false, xticklabelsize=lbsize, ygridvisible=false)
    axτ1 = Axis(fig[1, 1], 
        xaxisposition=:top, xlabel="Maintenance Date", xticks=(τ, [L"\tau_{%$i}" for i in eachindex(τ)]), xgridvisible=false, xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false)
    linkxaxes!(ax1, axτ1)
    linkyaxes!(ax1, axτ1)

    # New minimal x-axis separator
    ax_sep = Axis(fig[2, 1], 
        xlabel="Inspection Date", xticks=(inspection_dates[I], [L"t_{1, 1}", L"t_{1, 2}", L"t_{2, 1}", L"t_{4, 1}", L"t_{5, 1}", L"t_{5, 2}", L"t_{5, 3}"]), xticklabelsize=lbsize, xlabelsize=lbsize2,
        yticks=([], []), yticklabelsvisible=false, ygridvisible=false, xgridvisible=false,
        height=80)  # Small height for separator
    hidespines!(ax_sep, :l, :r, :t)
    ax_sep.yautolimitmargin = (0, 0)
    linkxaxes!(ax1, ax_sep)

    # Axis definition for the second figure (now in row 3)
    ax2 = Axis(fig[3, 1], xticks=(inspection_dates[I_char2], [L"t_{\nu_1^{(2)}, \phi_{1, 1}^{(2)}}", L"t_{\nu_1^{(2)}, \phi_{1, 2}^{(2)}}", L"t_{\nu_2^{(2)}, \phi_{2, 1}^{(2)}}", L"t_{\nu_2^{(2)}, \phi_{2, 2}^{(2)}}"]), ylabel = L"\textrm{Characteristic } (2): ARD_\infty", ylabelsize=lbsize2, xgridvisible=false, xticklabelsize=lbsize, ygridvisible=false)
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    # Set row height for separator
    rowsize!(fig.layout, 2, Relative(0.08))

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linewidth=1.)
    char2 = lines!(ax2, inspection_dates, Y∞, color=:red, linewidth=1.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax2, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax2, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax2, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax2, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 35.)

    # Add observation points with red color and bigger markersize
    observations_char1 = scatter!(ax1, inspection_dates[I_char1], Y1[I_char1], color=:black, markersize=20.)
    observations_char2 = scatter!(ax2, inspection_dates[I_char2], Y∞[I_char2], color=:black, markersize=20.)
    vlines!(ax1, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)
    vlines!(ax2, inspection_dates[I], color=:grey, linestyle=:dash, linewidth=1., alpha=0.5)

    # Legend
    axislegend(ax1, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(1)} = 0.5", L"\rho_2^{(1)} = 0.8", L"\rho_3^{(1)} = 0.9"], position = :lt, labelsize=lbsize2)
    axislegend(ax2, merge = true, [maint_type1, maint_type2, maint_type3], [L"\rho_1^{(2)} =0.2", L"\rho_2^{(2)} =0.5", L"\rho_3^{(2)} = 0.9"], position = :lt, labelsize=lbsize2)

    # virtual increments for the sole last observations of the first characteristic
    ps1=Point2f.([(t, 0) for t in vcat(0., inspection_dates[I_char1])])
    text_obs_var = [L"Z_0^{(1)}", L"Z_1^{(1)}", L"Z_2^{(1)}", L"\Delta Y_{3, 2}"]
    for i in 1:length(ps1)-1
        arrows2d!(ax1, ps1[i], ps1[i+1],
        argmode = :endpoint,
        tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8
        )
        text!(ax1, (ps1[i][1] + ps1[i+1][1]) / 2 - .5, 0., text=text_obs_var[i], color=:black, fontsize=lbsize - 5)
    end

    # virtual increments for the sole last observations of the second characteristic
    ps2=Point2f.([(t, 0) for t in vcat(0., inspection_dates[I_char2])])
    text_obs_var = [L"\Delta X_{1, 1}^{(2)}", L"\Delta Y_{1, 2}^{(2)}", L"Z_1^{(2)}", L"\Delta Y_{2, 2}^{(2)}"]
    for i in 1:length(ps2)-1
        arrows2d!(ax2, ps2[i], ps2[i+1],
        argmode = :endpoint,
        tail = Point2f[(0, 0), (1, -0.5), (1, 0.5)], taillength = 8
        )
        text!(ax2, (ps2[i][1] + ps2[i+1][1]) / 2 - .5, 0., text=text_obs_var[i], color=:black, fontsize=lbsize - 5)
    end

    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure_observed_increments_illustration.png", fig)


# One figure with both characteristics
fig = Figure(resolution = (1920, 1080))

    # Axis definition
    ax1 = Axis(fig[1, 1], ylabel = "Degradation level", ylabelsize=lbsize2, xlabel = "Time", xlabelsize=lbsize2, xgridvisible=false)

    # Main Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, X[1, :], color=:blue, linewidth=1.)
    char2 = lines!(ax1, inspection_dates, X[2, :], color=:red, linewidth=1.)
    # ARD Wiener process trajectory
    char1 = lines!(ax1, inspection_dates, Y1, color=:blue, linestyle=:dash, linewidth=1.)
    char2 = lines!(ax1, inspection_dates, Y∞, color=:red, linestyle=:dash, linewidth=1.)


    # Maintenance dates illustration with dashed colored vertical lines
    maint_type1 = vlines!(ax1, τ[1], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type2 = vlines!(ax1, τ[2], linestyle=:dash, color=:purple, linewidth=2.)
    vlines!(ax1, τ[3], linestyle=:dash, color=:magenta, linewidth=2.)
    maint_type3 = vlines!(ax1, τ[4], linestyle=:dash, color=:brown, linewidth=2.)

    ylims!(ax1, -2.5, 40.)

    vlines!(ax1, inspection_dates[I], color=:black, linewidth=1., linestyle=:dash)
I
    # Display and save the figure
    display(fig)
    save("C:\\Users\\bergerto\\Documents\\PhD\\Latex\\pics\\mw_main_figure_laurent_version.png", fig)