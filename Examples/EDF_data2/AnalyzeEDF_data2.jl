using CSV
using ARDmodels
using CairoMakie
using DataFrames

inch = 96
pt = 4/3
cm = inch / 2.54

systs = Vector{DegradationData}(undef, 9)

for i in 1:length(systs)
    maint = DataFrame(CSV.File(string("ARDmodels/exemples/EDF_data2/maint",i,".csv")))
    maint.TYPE = convert(Vector{String}, maint.TYPE)
    deg = DataFrame(CSV.File(string("ARDmodels/exemples/EDF_data2/deg",i,".csv")))
    infos = Dict("name" => string("Syst ", i))
    systs[i] = DegradationData(maint, deg,infos)
end

# 3 different maintenance types
maintTypes = sort(unique(reduce(vcat,[maintenances(systs[i]).TYPE for i in 1:length(systs)])))

f = Figure(size=(20cm, 20cm))
plot!(f, systs,
    maintenanceslinestyles = DataFrame(TYPE = maintTypes, LINESTYLE = [:solid, :dash, :dot]),
    markersize = 20, markercolor = :red, linecolor=:red, nbcols=3)
f
save("ARDmodels/exemples/EDF_DATA2/Datas.pdf",f)

#We define a WienerARD∞ model with the corresponding maintenance type
m = WienerARD∞(maintTypes)
# We compute the maximum likelihood estimator on the complet data set
fit_mle!(m,systs)
θ = params(m)
# And also for each signle trajectory
θs = DataFrame(Matrix{Union{Float64,Missing}}(missing, length(systs)+1, ncol(params(m))), :auto)
rename!(θs, vcat(names(params(m))))
θs[1, 1:ncol(θ)] = θ[1, :]
θs.syst .= "all"
for i in 1:length(systs)
    θs[i+1, 1:ncol(θ)] = fit_mle(m,systs[i])[1,:]
    θs.syst[i+1] = string(i)
end
print(θs)

#It is also possible to compute asymptotic confidence intervals instead of ponctual estimations
res = confint(m, systs)
ints = DataFrame(Matrix{Union{Float64, Missing}}(undef, length(systs) +1 , 2*ncol(params(m))), :auto)
ints.syst .= ""
rename!(ints,vcat(reduce(vcat, [string.(names(res)[i], ["_Inf", "_Sup"]) for i in 1:ncol(res)]), ["syst"]))
ints[1, 1:2:(2*ncol(params(m))-1)] = collect(res[1,:])
ints[1, 2:2:2*ncol(params(m))] = collect(res[2,:])
ints[1, 2*ncol(params(m))+1] = "all"
for i in 1:length(systs)
    res = confint(m, systs[i])
    ints[i+1, 1:2:(2*ncol(params(m))-1)] = collect(res[1,:])
    ints[i+1, 2:2:2*ncol(params(m))] = collect(res[2,:])
    ints[i+1, 2*ncol(params(m))+1] = string(i)
end
println(ints)

#We can simulate equivalent data sets (same number of systems, dates of maintenances and degradation observations)
#according to the parameters value estimate with the complet data set
#then, it is possible to evaluate the effective coverage rate of asympotic confidence intervals
#and also the quality of parameters estimations
#in the cases where estimation is done over all systems or only one of the system
datasim = deepcopy(systs)
θ = collect(params(m)[1,:])
nbMC = 500
θestMC = Array{Union{Float64, Missing}}(missing, nbMC, length(systs)+1, length(θ) )
CRMC = convert(Matrix{Union{Float64, Missing}}, fill(0., length(systs)+1, length(θ)))
for i in 1: nbMC
    print(string("=", i, "="))
    rand!(m, datasim)
    θestMC[i, 1, 1:length(θ)] =  collect(fit_mle(m, datasim)[1,:])
    res = confint(m, datasim)
    for j in 1:length(θ)
        CRMC[1, j] += (res[1, j] <= θ[j] <= res[2, j]) ? 1 : 0
    end
    for k in 1:length(systs)
        θestMC[i, k+1, 1:length(θ)] =  collect(fit_mle(m, datasim[k])[1,:])
        res = confint(m, datasim[k], print_summary=false)
        for j in 1:length(θ)
            if ismissing(res[1, j])
                CRMC[k+1, j] = missing
            else
                CRMC[k+1, j] += (res[1, j] <= θ[j] <= res[2, j]) ? 1 : 0
            end 
        end
    end
end
f = Figure(size=(80cm, 80cm))
axs = Array{Axis}(undef, length(θ), length(systs)+1)
CR = DataFrame(Matrix{Union{Float64, Missing}}(missing, length(systs)+1, length(θ)), :auto)
rename!(CR, names(params(m)))
CR.syst .= ["all"]
xmins = [-10e-5, 0., -0.5, -1., -1.]
xmaxs = [60e-5, 0.08, 1.5, 2., 1.5]
for j in 1:length(θ)
    CR[1, j] = CRMC[1, j] / nbMC
    axs[j, 1] = Axis(f[j, 1], title= string(names(params(m))[j]," for all"))
    Makie.hist!(axs[j, 1], θestMC[:,1 ,j], normalization = :pdf)
    Makie.vlines!(axs[j, 1], θ[j], color=:black)
    Makie.xlims!(axs[j, 1], xmins[j], xmaxs[j])
end
for i in 1:length(systs)
    CR[i+1, 1:length(θ)] = CRMC[i+1, 1:length(θ)] / nbMC
    CR.syst[i+1] = string(i)
    for j in 1:length(θ)
        axs[j, i+1] = Axis(f[j, i+1], title= string(names(params(m))[j]," for syst ", i))
        if !ismissing(θestMC[1, i+1 ,j])
            Makie.hist!(axs[j, i+1], θestMC[:, i+1 ,j], normalization = :pdf)
        end
        Makie.vlines!(axs[j, i+1], θ[j], color=:black)
        Makie.xlims!(axs[j, i+1], xmins[j], xmaxs[j])
    end
end
f
save("ARDmodels/exemples/EDF_DATA2/DistribParamEstim.pdf",f)
println(CR)
#10×6 DataFrame
# Row │ μ            σ            ρ_P0         ρ_P1         ρ_P2         syst   
#     │ Float64?     Float64?     Float64?     Float64?     Float64?     String 
#─────┼─────────────────────────────────────────────────────────────────────────
#   1 │       0.912        0.888        0.952        0.908        0.908  all
#   2 │       0.84         0.64         0.742        0.778        0.792  1
#   3 │ missing      missing      missing      missing      missing      2
#   4 │       0.744        0.554        0.742        0.78   missing      3
#   5 │       0.786        0.644  missing      missing            0.81   4
#   6 │       0.618        0.448        0.722        0.676  missing      5
#   7 │       0.746        0.588        0.75   missing            0.782  6
#   8 │       0.692        0.608        0.94   missing            0.768  7
#   9 │ missing      missing      missing      missing      missing      8
#  10 │       0.37         0.274        0.6          0.472        0.494  9

#Finally, we can easily plot the quantiles of the distributions of the degradation values
#depending on the fact that 
# - we only know that degradation is equal to 0 at time 0
# - or given the last degradation observation

f = Figure(size=(20cm,20cm))
plot!(f, systs,
    maintenanceslinestyles = DataFrame(TYPE = maintTypes, LINESTYLE = [:solid, :dash, :dot]),
    markersize = 20, markercolor = :red, linecolor=:red, nbcols=3)
for i in 1:length(systs)
    plot!(f.content[i ,1],m, systs[i], QuantilesPlot, tmaxs=[3e5], color=:blue)
    plot!(f.content[i ,1],m , systs[i], QuantilesPlot, froms=convert(Vector{Int64}, 0:nrow(systs[i].degradations)), tmaxs=vcat(systs[i].degradations.DATE, 3e5), color=:red)
    Makie.ylims!(f.content[i ,1],0, 40)
end
f
save("ARDmodels/exemples/EDF_DATA2/DegObsQuantiles.pdf",f)

#We can forcast the futur of the degradation process given each of the last observed degradation values of the degradation
# and in particular also forcast the effect of futur maintenances
Syst8Frocast = deepcopy(systs[8])
push!(Syst8Frocast.maintenances, [2.8e5, maintTypes[2]])#suppose that for system 8 a maintenance of type P1 is planned at time 3.e5
f=Figure()
plot!(f, Syst8Frocast, 
    maintenanceslinestyles=DataFrame(TYPE = maintTypes, LINESTYLE = [:solid, :dash, :dot]),
    markercolor=:red, linecolor=:red)
f#look at the new planned maintenances after the latest degradation observations
plot!(f.content[1], m , Syst8Frocast, QuantilesPlot, froms=[3], tmaxs=[5.e5], color=:red)
plot!(f.content[1], m , Syst8Frocast, QuantilesPlot, froms=[2], tmaxs=[5.e5], color=:blue)
plot!(f.content[1], m , Syst8Frocast, QuantilesPlot, froms=[1], tmaxs=[5.e5], color=:green)
Makie.ylims!(f.content[1], 0, 5.e5)
Makie.ylims!(f.content[1], 0, 60.)
f#we are able to forcast the condtionnal futur evolution of the degradation using the estimated model

