using DataFrames
using ARDmodels
using CSV
using CairoMakie

systs = Vector{DegradationData}(undef, 9)


for i in 1:9
    maint = DataFrame(CSV.File(string("/Users/laurentdoyen/OwnCloud-2/Recherches/JuliaProg/ARDmodels/exemples/dataExemple/maint",i,".csv")))
    maint.TYPE = convert(Vector{String}, maint.TYPE)
    deg = DataFrame(CSV.File(string("/Users/laurentdoyen/OwnCloud-2/Recherches/JuliaProg/ARDmodels/exemples/dataExemple/deg",i,".csv")))
    infos = Dict("name" => string("Syst ", i))
    systs[i] = DegradationData(maint, deg,infos)
end

println(systs)

inch = 96
pt = 4/3
cm = inch / 2.54

CairoMakie.activate!()

f = Figure()
ax = Axis(f[1,1])
plot!(ax, systs[1], MaintenancesPlot)
f
f = Figure()
ax = Axis(f[1,1])
plot!(ax, systs[1], DegradationsPlot)
f
f = Figure()
ax = Axis(f[1,1])
plot!(ax, systs[1])
f
f = Figure()
ax = Axis(f[1,1])
plot!(ax, systs[1],maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]))
f
f = Figure()
ax = Axis(f[1,1])
plot!(ax, systs[1],
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]),
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot])
    )
f
systs[1]
f = Figure()
plot!(f, systs[1],
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]),
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot])
    )
f
f = Figure()
plot!(f, systs, 
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]), 
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot]),
    nbcols=3)
f

systs[1].degradations.CENS[3]=1
f = Figure()
plot!(f, systs, 
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]), 
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot]),
    nbcols=3,
    plot_cens=true, plot_pb_nb_maint=true,
    markercolor=:red, markersize=10)
f
systs[1].degradations.CENS[3]=0

f = Figure()
plot!(f, systs, 
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]), 
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot]),
    nbcols=3,
    plot_cens=true, plot_pb_nb_maint=true,
    marker='x',
    linecolor=:red, linestyle=:solid,
    xlims=[-1. *10^5, 3. *10^5], ylims=[-10.,40.])
f

f = Figure()
ax = Axis(f[1,1])
plot!(ax, systs[1], linestyle=:none,
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]),
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot])
)
f
m = WienerARD∞(["C", "P1", "P2"], 2.4e-4, 0.031, [0.8, 0.16, -0.1])
plot!(ax, m, systs[1], QuantilesPlot, step = 10)
f
plot!(ax, m, systs[1], QuantilesPlot, [0.5], tmaxs=[2.5e5], linestyles=[:solid])
f
plot!(ax, m, systs[1], QuantilesPlot, froms=[5], color=:green)
f
plot!(ax, m, systs[1], QuantilesPlot, tmaxs=[2.5e5], froms=[4], color=:blue)
f
plot!(ax, m, systs[1], QuantilesPlot, [0.5], tmaxs=[2.5e5],froms=[4], linestyles=[:solid], color=:blue)
f

syst = deepcopy(systs[1])
syst
syst.degradations = syst.degradations[[1, 2, 3, 4, 5, 7, 8, 9, 10, 12], :]
f = Figure()
ax = Axis(f[1,1])
plot!(ax, syst, linestyle=:none,
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]),
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot])
)
f
plot!(ax, m, syst, QuantilesPlot, tmaxs=[syst.degradations.DATE[8]])
f
plot!(ax, m, syst, QuantilesPlot, froms=[6], color=:green)
f
plot!(ax, m, syst, QuantilesPlot, froms=[4], tmaxs=[syst.degradations.DATE[8]], color=:blue)
f

f = Figure()
ax = Axis(f[1,1])
plot!(ax, systs[1], linestyle=:none,
    maintenancescolors = DataFrame(TYPE =["C", "P1", "P2"], COLOR = [:red, :green, :blue]),
    maintenanceslinestyles = DataFrame(TYPE =["C", "P1", "P2"], LINESTYLE = [:solid, :dash, :dot])
)
f
plot!(ax, m, systs[1], QuantilesPlot, 
    froms=[0], 
    tmaxs=[2.5e5])
f
plot!(ax, m, systs[1], QuantilesPlot, 
    froms=convert(Vector{Int64}, 1:nrow(systs[1].degradations)), 
    tmaxs=vcat(systs[1].degradations.DATE[2:nrow(systs[1].degradations)], 2.5e5 ), color=:green)
f
plot!(ax, m, systs[1], QuantilesPlot, 
    froms=[9], 
    tmaxs=[2.5e5 ], color=:blue)
f