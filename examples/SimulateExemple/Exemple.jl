using ArithmeticReductionDegradationModels
using CairoMakie
using DataFrames
using Random

inch = 96
pt = 4/3
cm = inch / 2.54

step = 0.01 #simulation step

#Comparison ARD1 and ARD∞
m1 = WienerARD1(["P"], 2., 1.5, [0.5])#maintenance reduces degradation of on half of the suppelment of degradation accumulated since the last maintenance
m∞ = WienerARD∞(["P"], 2., 1.5, [0.5])#maintenance reduces degradation of an half of its value just before maintenance
degs = Vector{DegradationData}(undef, 2)
Random.seed!(123)
degs[1] = rand(m1, collect(1.:10.), step, 11., fill("P", 10))
infos!(degs[1], Dict("name" => "ARD1"))
Random.seed!(123)#in this case both ARD1 and ARD∞ are constructs from the same underlying degradation
degs[2] = rand(m∞, collect(1.:10.), step, 11., fill("P", 10))
infos!(degs[2], Dict("name" => "ARD∞"))
f=Figure()
plot!(f, degs, 
    marker=:none, linestyle=:solid, linecolor=(:blue, 0.5))
f

#ARD(∞)
mT = ["P1", "P2"]#two different maintenances
m = WienerARD∞(mT, 2., 1.5, [0.4, 0.8])#maintenance P2 is twice more efficient than P1

#2 independent systems following the same model but with different maintenances
maints = Vector{DataFrame}(undef, 2)
maints[1] = DataFrame(
    DATE = [1., 2.5, 3.5],
    TYPE = [mT[1], mT[1], mT[2]])
maints[2] = DataFrame(
    DATE = [2.],
    TYPE = [mT[2]])
endTimes = [4.5, 3]

degs = Vector{DegradationData}(undef, 2)
for i in 1:length(degs)
    degs[i] = rand(m, maints[i].DATE, step, endTimes[i], maints[i].TYPE)
    infos!(degs[i], Dict("name" => string("ARD∞ syst", i)))
end
f=Figure()
plot!(f, degs,
    maintenanceslinestyles=DataFrame(TYPE = mT, LINESTYLE = [:solid, :dot]),
    marker=:none, linestyle=:solid, linecolor=RGBf(1., 0., 0.))
f#One simulation of the evolution of the degradation
nbsim = 20
for i in 2:(nbsim)
    rand!(m, degs)
    for j in 1: length(degs)
        plot!(f.content[j], degs[j],
            maintenanceslinestyles=DataFrame(TYPE = mT, LINESTYLE = [:none, :none]),
            marker=:none, linestyle=:solid, linecolor=RGBf((nbsim-i+1)/(nbsim-1), 0., (i-1.)/(nbsim-1)))
    end
end
for j in 1: length(degs)
    Makie.ylims!(f.content[j],-2, 8)
end
f#20 iid simulation of the evolution of the degradation
for i in 1:length(degs)
    plot!(f.content[i], m , degs[i], QuantilesPlot, tmaxs=[endTimes[i]+2.], color=:green)
end
f#at each time degradation is normally ditributed, the green band theoritcally contains 50% of the observed degradation values

f=Figure()
plot!(f, degs,
    maintenanceslinestyles=DataFrame(TYPE = mT, LINESTYLE = [:solid, :dot]),
    marker=:none, linestyle=:solid, linecolor=RGBf(1., 0., 0.))
f#Strating from one simulated degradation
from = 150#We simulate the conditional evolution of the degradation given its vlaue at time from*step=150*0.01.1.5
nbsim = 20
for i in 2:(nbsim)
    for j in 1: length(degs)
        rand!(m, degs[j], from)
        plot!(f.content[j], degs[j],
            maintenanceslinestyles=DataFrame(TYPE = mT, LINESTYLE = [:none, :none]),
            marker=:none, linestyle=:solid, linecolor=RGBf((nbsim-i+1)/(nbsim-1), 0., (i-1.)/(nbsim-1)))
    end
end
for j in 1: length(degs)
    Makie.ylims!(f.content[j],-2, 8)
end
f
for i in 1:length(degs)
    plot!(f.content[i], m , degs[i], QuantilesPlot, froms=[from], tmaxs=[endTimes[i]+2.], color=:green)
end
f#the condition value of the degradation is still normally distributed and we are still able to compute its caracteristics: ear a 50% band

#Degradation observations and inference
#lets take a trajectory
f=Figure()
plot!(f, degs,
    maintenanceslinestyles=DataFrame(TYPE = mT, LINESTYLE = [:solid, :dot]),
    marker=:none, linestyle=:solid, linecolor=(:blue, 0.3))
f
#We suppose that we only know the degradation value at some specific times
#println(degs∞[1].degradations)
indObs1 = [101, 138, 252, 302, 356, 455]
#println(degs∞[2].degradations)
indObs2 = [142, 190, 216, 301]
datas = Vector{DegradationData}(undef, 2)
datas[1] = DegradationData(
    maintenances(degs[1]), 
    degradations(degs[1])[indObs1, :],
    Dict("name"=>"syst1")) 
datas[2] = DegradationData(
        maintenances(degs[2]), 
        degradations(degs[2])[indObs2, :],
        Dict("name"=>"syst2")) 
for i in 1:length(datas)
    plot!(f.content[i], datas[i], markercolor=:red, linecolor=:red)
end
f

mest = deepcopy(m)
fit_mle!(mest, datas)
θest = params(mest)#The estimated parameters value knowing only the degradation at times plotted in red
θest#The estimated parameters value knowing only the degradation at times plotted in red
θest .- params(m)#The biais
confint(mest, datas)#asympotic confidence intervals for parameter values

f=Figure()
plot!(f, vcat(degs, degs), 
    maintenanceslinestyles=DataFrame(TYPE = mT, LINESTYLE = [:solid, :dot]),
    marker=:none, linestyle=:solid, linecolor=(:blue, 0.5), nbcols=2)
f
for i in 1:length(datas)
    plot!(f.content[i], datas[i], markercolor=:red, linecolor=:red)
    plot!(f.content[i], m , datas[i], QuantilesPlot, tmaxs=[endTimes[i]+2.], color=:blue)
    plot!(f.content[i], m , datas[i], QuantilesPlot, froms=convert(Vector{Int64}, 0:nrow(datas[i].degradations)), tmaxs=vcat(datas[i].degradations.DATE, endTimes[i]+2.), color=:red)
    Makie.ylims!(f.content[i],0, 6)
    plot!(f.content[i+2], datas[i], markercolor=:red, linecolor=:red)
    plot!(f.content[i+2], mest , datas[i], QuantilesPlot, tmaxs=[endTimes[i]+2.], color=:blue)
    plot!(f.content[i+2], mest , datas[i], QuantilesPlot, froms=convert(Vector{Int64}, 0:nrow(datas[i].degradations)), tmaxs=vcat(datas[i].degradations.DATE, endTimes[i]+2.), color=:red)
    Makie.ylims!(f.content[i+2],0, 6)
end
f#Comparison of 50% confidence bands and conditional confidence bands using the real parameter values or the estimated ones

#Forecasting the futur evolution of the degradation
DataFrocast = deepcopy(datas)
push!(DataFrocast[1].maintenances, [5., mT[2]])#suppose that for system 1 a maintenance of type P2 is planned at time 5.
push!(DataFrocast[2].maintenances, [4., mT[1]])
push!(DataFrocast[2].maintenances, [5., mT[2]])#and system 2, two maintenances at times 4. and 5. of types P1 and P2 are planned
f=Figure()
plot!(f, DataFrocast, 
    maintenanceslinestyles=DataFrame(TYPE = mT, LINESTYLE = [:solid, :dot]),
    markercolor=:red, linecolor=:red)
f#look at the new planned maintenances after the latest degradation observations
plot!(f.content[1], mest , DataFrocast[1], QuantilesPlot, froms=[6], tmaxs=[6.], color=:red)
plot!(f.content[2], mest , DataFrocast[2], QuantilesPlot, froms=[4], tmaxs=[6.], color=:red)
Makie.xlims!(f.content[1], 0, 6.)
Makie.xlims!(f.content[2], 0, 6.)
f#we are able to forcast the condtionnal futur evolution of the degradation using the estimated model
plot!(f.content[1], m , DataFrocast[1], QuantilesPlot, froms=[6], tmaxs=[6.], color=:blue)
plot!(f.content[2], m , DataFrocast[2], QuantilesPlot, froms=[4], tmaxs=[6.], color=:blue)
f# since this simulated data it can compared with the real one using real parameter values instead of estimated ones

#Finally the plot values can also be explicitly obtained 
QData2Frocast = deepcopy(DataFrocast[2])
QData2Frocast
push!(QData2Frocast.degradations, [6., 3, 0., ""])
forcastQ = quantile(mest , QData2Frocast, [0.25, 0.5, 0.75], from=4)[5, :]

