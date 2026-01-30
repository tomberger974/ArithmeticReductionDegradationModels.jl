using DataFrames
using CairoMakie
using ArithmeticReductionDegradationModels

m = WienerARD∞()
params(m)
θ = params(m)
θ[1,2] = 1.3
θ 
params(m)
m = WienerARD∞("")
params(m)
m = WienerARD∞("m")
params(m)
m = WienerARD∞(["m1"])
params(m)
m = WienerARD∞(["m1", "m2"])
params(m)
m = WienerARD1()
params(m)
m = WienerARD1("")
params(m)
m = WienerARD1("m")
params(m)
m = WienerARD1(["m1"])
params(m)
m = WienerARD1(["m1", "m2"])
params(m)



m = Vector{WienerARD∞}(undef, 2)

i = 1
m[i] = WienerARD∞("Maint", 2., 1.5, 0.8)
params(m[i])
params!(m[i], [2.3, 1.57,0.88])
params(m[i])
params!(m[i], [2.3, 1.57,0.88], ["maint"])
params(m[i])

i = 2
m[i] = WienerARD∞(["M1", "M2"], 2., 1.5, [0.5, 0.8])
params(m[i])
params!(m[i], [2.3, 1.57, 0.55, 0.86])
params(m[i])
params!(m[i], [2.3, 1.57,0.54, 0.89], ["m1", "m2"])
params(m[i])

minf = deepcopy(m)
m = Vector{WienerARD1}(undef, 2)

i = 1
m[i] = WienerARD1("Maint", 2., 1.5, 0.8)
params(m[i])
params!(m[i], [2.3, 1.57,0.88])
params(m[i])
params!(m[i], [2.3, 1.57,0.88], ["maint"])
params(m[i])
collect(params(m[i])[1,:])

i = 2
m[i] = WienerARD1(["M1", "M2"], 2., 1.5, [0.5, 0.8])
params(m[i])
params!(m[i], [2.3, 1.57, 0.55, 0.86])
params(m[i])
params!(m[i], [2.3, 1.57,0.54, 0.89], ["m1", "m2"])
params(m[i])

typeof(m[2]) == WienerARD1

deg = DataFrame(DATE = [0.6, 1.2, 1.4, 1.3, 2.5, 2.2, 4.0, 5.0, 5.0, 5.2], 
    NB_MAINTENANCES = [0, 1, 1, 1, 2, 2, 4, 4, 5, 5],
    VALUE = fill(3.1, 10))
maint = DataFrame(DATE = [1., 2., 3., 4., 5.])
data = DegradationData(maint, deg)
degradations(data)
maintenances(data)
d = degradations(data)
d.DATE[2] = 1.3
d
degradations(data)
maint = maint[[3,2,5,1,4],:]
maintenances!(data, maint)
println(maintenances(data))
deg = deg[[3, 4, 9, 7, 5, 1, 2, 10, 6, 8], :]
degradations!(data, deg)
println(degradations(data))
data = DegradationData(maint, deg)
println(degradations(data))
println(maintenances(data))

using Random
Random.seed!(1234)
m = WienerARD∞("", 2., 1.5, 0.8)
rand!(m,data)
println(degradations(data))
degradations(data).VALUE[8] * (1-0.8)
degradations(data).VALUE[9]
d = rand(m, 1., 5, 0.5)
println(degradations(d))
println(maintenances(d))

Random.seed!(1234)
m = WienerARD1("", 2., 1.5, 0.8)
rand!(m,data)
println(degradations(data))
degradations(data).VALUE[8] - 0.8 * (degradations(data).VALUE[8]-degradations(data).VALUE[7])
degradations(data).VALUE[9]

maint = DataFrame(DATE = [1., 2., 3., 4., 5.], TYPE = ["m1", "m1", "m2", "m1", "m2"])
m = WienerARD1(["m1","m2"], 2., 1.5, [0.8, 0.2])
maintenances!(data, maint)
rand!(m,data)
println(degradations(data))
degradations(data).VALUE[8] - 0.2 * (degradations(data).VALUE[8]-degradations(data).VALUE[7])
degradations(data).VALUE[9]


maint = DataFrame(DATE = [1., 2., 3., 4., 5.], TYPE = ["m1", "m1", "m2", "m1", "m2"])
m = WienerARD∞(["m1","m2"], 2., 1.5, [0.8, 0.2])
maintenances!(data, maint)
rand!(m,data)
println(degradations(data))
degradations(data).VALUE[8] *(1 - 0.2)
degradations(data).VALUE[9]

rand!(m,data,4)
println(degradations(data))