
using CSV
using ArithmeticReductionDegradationModels
using CairoMakie
using DataFrames

inch = 96
pt = 4/3
cm = inch / 2.54

systs = Vector{DegradationData}(undef, 4)

for i in 1:length(systs)
    maint = DataFrame(CSV.File(string("ARDmodels/exemples/EDF_data1/maint",i,".csv")))
    maint.TYPE = convert(Vector{String}, maint.TYPE)
    deg = DataFrame(CSV.File(string("ARDmodels/exemples/EDF_data1/deg",i,".csv")))
    infos = Dict("name" => string("Syst ", i))
    systs[i] = DegradationData(maint, deg,infos)
end

# Only one maintenance types
maintTypes = unique(reduce(vcat,[maintenances(systs[i]).TYPE for i in 1:length(systs)]))

f = Figure()
plot!(f, systs, markersize = 20, markercolor = :red, linecolor=:red, nbcols=2)
f
save("ARDmodels/exemples/EDF_DATA1/Datas.pdf",f)


#We define a WienerARD∞ model with the corresponding maintenance type
m = WienerARD∞(maintTypes)
# We compute the maximum likelihood estimator on the complet data set
fit_mle!(m,systs)
θ = params(m)
# And also for each signle trajectory
θs = deepcopy(θ)
θs.syst = ["all"]
for i in 1:length(systs)
    res = fit_mle(m,systs[i])
    res.syst = [string(i)]
    push!(θs, res[1,:])
end
print(θs)

#It is also possible to compute asymptotic confidence intervals instead of ponctual estimations
res = confint(m, systs)
ints = DataFrame(Matrix{Float64}(undef, length(systs) +1 , 2*ncol(params(m))), :auto)
ints.syst .= ""
rename!(ints,vcat(reduce(vcat, [string.(names(res)[i], ["Inf", "Sup"]) for i in 1:ncol(res)]), ["syst"]))
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
θestMC = Array{Float64}(undef, nbMC, length(systs)+1, length(θ) )
CRMC = fill(0., length(systs)+1, length(θ))
for i in 1: nbMC
    rand!(m, datasim)
    θestMC[i, 1, 1:length(θ)] =  collect(fit_mle(m, datasim)[1,:])
    res = confint(m, datasim)
    for j in 1:length(θ)
        CRMC[1, j] += (res[1, j] <= θ[j] <= res[2, j]) ? 1 : 0
    end
    for k in 1:length(systs)
        θestMC[i, k+1, 1:length(θ)] =  collect(fit_mle(m, datasim[k])[1,:])
        res = confint(m, datasim[k])
        for j in 1:length(θ)
            CRMC[k+1, j] += (res[1, j] <= θ[j] <= res[2, j]) ? 1 : 0
        end
    end
end
f = Figure(size=(30cm, 30cm))
axs = Array{Axis}(undef, length(θ), length(systs)+1)
CR = DataFrame(Matrix{Float64}(undef, length(systs)+1, length(θ)), :auto)
rename!(CR, names(params(m)))
CR.syst .= ["all"]
xmins = [2e-5, 0, 0]
xmaxs = [8e-5, 0.008, 0.5]
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
        Makie.hist!(axs[j, i+1], θestMC[:, i+1 ,j], normalization = :pdf)
        Makie.vlines!(axs[j, i+1], θ[j], color=:black)
        Makie.xlims!(axs[j, i+1], xmins[j], xmaxs[j])
    end
end
f
save("ARDmodels/exemples/EDF_DATA1/DistribParamEstim.pdf",f)
println(CR)

#Finally, we can easily plot the quantiles of the distributions of the degradation values
#depending on the fact that 
# - we only know that degradation is equal to 0 at time 0
# - or given the last degradation observation
f = Figure(size=(20cm,20cm))
plot!(f, systs, markersize = 15, markercolor = :red, linecolor=:red, nbcols=2)
for i in 1:length(systs)
    plot!(f.content[i ,1],m, systs[i], QuantilesPlot, tmaxs=[3e5], color=:blue)
    plot!(f.content[i ,1],m , systs[i], QuantilesPlot, froms=convert(Vector{Int64}, 0:nrow(systs[i].degradations)), tmaxs=vcat(systs[i].degradations.DATE, 3e5), color=:red)
    Makie.ylims!(f.content[i ,1],0, 14)
end
f
save("ARDmodels/exemples/EDF_DATA1/DegObsQuantiles.pdf",f)

