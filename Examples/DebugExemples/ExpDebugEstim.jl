using DataFrames
using ARDmodels
using CSV
using Distributions

systs = Vector{DegradationData}(undef, 9)

for i in 1:9
    maint = DataFrame(CSV.File(string("/Users/laurentdoyen/OwnCloud-2/Recherches/JuliaProg/ARDmodels/exemples/dataExemple/maint",i,".csv")))
    maint.TYPE = convert(Vector{String}, maint.TYPE)
    deg = DataFrame(CSV.File(string("/Users/laurentdoyen/OwnCloud-2/Recherches/JuliaProg/ARDmodels/exemples/dataExemple/deg",i,".csv")))
    infos = Dict("name" => string("Syst ", i))
    systs[i] = DegradationData(maint, deg,infos)
end
systs

m = WienerARD∞(["C", "P1", "P2"])
fit_mle(m, systs)
fit_mle(m, systs, ρ0s=[0.1, 0.1, 0.1])
θest = fit_mle(m, systs)
θest_syst = deepcopy(θest)
θest_syst.syst = ["all"]
#MLE système par système
for i in 1:length(systs)
    res = fit_mle(m, [systs[i]])
    res.syst = [string(i)]
    θest_syst = outerjoin(θest_syst, res, on=names(θest_syst) ∩ names(res), matchmissing=:equal)
end
println(θest_syst)
mm = deepcopy(m)
fit_mle!(mm, systs)

loglikelihood(m, systs, [0.002, 0.06, 0.8, 0.1, 0.2])
loglikelihood(m, systs, [0.002, 0.06, 0.8, 0.1, 0.2], withDerivatives = true)
α = 0.05 # niveau de confiance 1-α
u = quantile(Normal(),1-α/2)
In = loglikelihood(m, systs, collect(θest[1,:]), withDerivatives=true)[3]
sqrt_inv_In = sqrt(inv(-In))
intConfAss = deepcopy(θest)
intConfAss[1,:] = collect(θest[1,:]) .- u .* [sqrt_inv_In[i,i] for i in 1:size(sqrt_inv_In)[1]]
push!(intConfAss, collect(θest[1,:]) .+ u .* [sqrt_inv_In[i,i] for i in 1:size(sqrt_inv_In)[1]])
intConfAss
confint(m, systs)
confint(m, systs, 0.95)

println(ARDmodels.mean_var(mm, systs[1]))
println(ARDmodels.mean_var(mm, systs[1], from=4,to=10))
println(mean(mm,systs[1], from=4,to=10))
println(mean(mm,systs[1]))
println(var(mm,systs[1], from=4,to=10))
println(var(mm,systs[1]))
println(std(mm,systs[1], from=4,to=10))
println(std(mm,systs[1]))
println(quantile(mm,systs[1], from=4,to=10))
println(quantile(mm,systs[1]))
println(quantile(mm,systs[1], [0.25,0.75], from=4,to=10))
println(quantile(mm,systs[1], [0.25,0.75]))

moyenne = mean(mm,systs[1], from=4,to=10)
moyenne.DATE[2] = 3.5
moyenne
systs[1]

