using CairoMakie
using Distributions

using StatsBase
using Distributions
using Random

using Pkg
Pkg.activate(".")

using Revise
using ArithmeticReductionDegradationModels

include("../nearest_time.jl")

abstract type Process end

mutable struct Wiener <: Process
    drift::Float64
    dispertion::Float64

    function Wiener(drift::Float64 = 0, dispersion::Float64 = 1)
        new(drift, dispersion)
    end

    function Wiener(underlyingDegradation::Vector{Float64})
        if length(underlyingDegradation) != 2
            error("The vector of underlying degradation contains only 2 parameteres")
        end
        new(underlyingDegradation[1], underlyingDegradation[2])
    end
end

#####TRAJECTOIRE D'UN MOUVEMENT BROWNIEN#####
function Brownian_trajectory(k=100, T=1)
    #creation of a time table plus a sample of the distribution 
    time = convert(Vector{Float64}, range(0, T, k))
    X = rand(Normal(0, 1), k)

    #simulation of the Brownian motion
    B = zeros(k)
    B[2] = sqrt(time[2])*X[2]#set B[1]=0 and initialization at B[2]
    for i ∈ 3:k
        B[i] = B[i-1] + sqrt(time[i] - time[i-1])*X[i]
    end

    return time, B
end

#####TRAJECTOIRE D'UN PROCESSUS DE WIENER#####
function trajectory(wiener::Wiener=Wiener(), k=100, T=1)
    μ, σ = wiener.drift, wiener.dispertion
    W = zeros(k)

    #On utilise un mouvement brownien
    time, B = Brownian_trajectory(k, T)
    
    #On en déduit le processus de Wiener
    W = μ*time .+ σ*B
    
    return time, W
end

#####TRAJECTOIRE D'UN PROCESSUS AUQUEL ON A AJOUTE UN UNIQUE TEMPS DE MAINTENANCE ARD∞#####
function trajectory(wienerARD∞::WienerARD∞=WienerARD∞(), τ::Float64=.5, k=100, T=1)
    if !(0 <= τ <= T)
        error("τ must be between 0 and T")
    end

    #simulate first a trajectory for the underlying process
    wiener = Wiener(wienerARD∞.underlyingDegradation)
    time, Y = trajectory(wiener, k, T)

    #approximate τ by the the first index of time above τ
    time_stop_index = 1
    while τ > time[time_stop_index]
        time_stop_index += 1
    end

    #after the time stop add the ARD∞ condition
    ρ = wienerARD∞.maintenances[1 ,2]
    for i in time_stop_index+1:k
        Y[i] = Y[i] - ρ*Y[time_stop_index]
    end

    return time, Y
end

function trajectory(wienerARD1::WienerARD1=WienerARD1(), τ::Float64=.5, k=100, T=1)
    if !(0 <= τ <= T)
        error("τ must be between 0 and T")
    end

    #simulate first a trajectory for the underlying process
    wiener = Wiener(wienerARD1.underlyingDegradation)
    time, Y = trajectory(wiener, k, T)

    #approximate τ by the the first index of time above τ
    time_stop_index = nearest_time(τ, time)

    #after the time stop add the ARD∞ condition
    ρ = wienerARD1.maintenances[1 ,2]
    for i in time_stop_index+1:k
        Y[i] = Y[i] - ρ*Y[time_stop_index]
    end

    return time, Y
end

function trajectory(maintenances::DataFrame, wienerARD::WienerARD, k=100, T=1)
    τ = maintenances.DATE

    for t ∈ τ 
        if !(0 <= t <= T)
            error("all elements of τ must be between 0 and T")
        end
    end

    #simulate first a trajectory for the underlying process
    wiener = Wiener(wienerARD.underlyingDegradation)
    time, Y = trajectory(wiener, k, T)

    #approximate τ by the the first index of time above τ[i]
    time_stop_indices = nearest_time(τ, time)

    #after the time stop add the ARD1 condition
    type_to_ρ = Dict(wienerARD.maintenances.TYPE .=> wienerARD.maintenances.VALUE)

    if typeof(wienerARD) == WienerARD1
        for i in eachindex(time_stop_indices)
            if i == 1
                a = Y[time_stop_indices[i]]
                ρ = type_to_ρ[maintenances[i, "TYPE"]]
                Y[time_stop_indices[i]:end] .-= ρ*a
            else
                a = Y[time_stop_indices[i]] - Y[time_stop_indices[i-1]]
                ρ = type_to_ρ[maintenances[i, "TYPE"]]
                Y[time_stop_indices[i]:end] .-= ρ*a
            end
        end
    else
        for i in eachindex(time_stop_indices)
            a = Y[time_stop_indices[i]]
            ρ = type_to_ρ[maintenances[i, "TYPE"]]
            Y[time_stop_indices[i]:end] .-= ρ*a
        end
    end

    return time, Y
end



#####GENERATEUR DE DATASET#####
function CreateData2(maintenanceDates::Vector{Float64}, Δt::Float64, t_max::Float64=maintenanceDates[end], maintenanceTypes::Union{Vector{String}, Vector{Symbol}}=fill("",length(maintenanceDates)); imaint_from::Int64=0, t_from::Float64=0.)

    #to avoid errors
    if Δt <= 0.
        error("Time increment Δt must be positive")
    end
    if length(maintenanceDates) != length(maintenanceTypes)
        error("Vectors of maintenancesDates and maintenancesTypes must have the same length.")
    end

    #initialize vectors
    degDates = Vector{Float64}()
    degNbMaint = Vector{Int64}()
    degTypes = Vector{String}()

    #In case where we want our observations to begin at maintenance imaint_from
    if imaint_from == 0
        maint = vcat(0,maintenanceDates[1:end])
    else
        maint = maintenanceDates[imaint_from:end]
    end

    for (jj,τ) in enumerate(maint)
        
        #we don't want to go farther than t_max
        if τ > t_max 
            break
        end
        
        #reindex j so that we begin at i_maint_from
        j = jj + imaint_from

        #create the "After" rows
        if (τ >= t_from)
            degDates = vcat(degDates, τ)
            degNbMaint = vcat(degNbMaint, j-1)
            degTypes = vcat(degTypes, "After")
        end


        #original code
        #t = max(t_from, τ + Δt)

        #my code
        if jj == 0
            t = max(t_from, τ + Δt)
        else
            t = t_from
            while t <= τ
                t += Δt
            end
        end


        τ_suiv = j <= length(maintenanceDates) ? maintenanceDates[j] : t_max
        while t < min(τ_suiv, t_max)
            degDates = vcat(degDates, t)
            degNbMaint = vcat(degNbMaint, j-1)
            degTypes = vcat(degTypes, "Between")
            t += Δt
        end
        if j <= length(maintenanceDates)
            if (τ_suiv >= t_from) && (τ_suiv <= t_max)
                if maintenanceDates[j] < t_max
                    degTypes = vcat(degTypes,  "Before" )
                else
                    degTypes = vcat(degTypes,  "Between" )
                end
                degDates = vcat(degDates, τ_suiv)
                degNbMaint = vcat(degNbMaint, j-1)
            end
        else
            if (maintenanceDates[end] < t_max) && (t >= t_from) && (t <= t_max)
                degTypes = vcat(degTypes,  "Between" )
                degDates = vcat(degDates, t)
                degNbMaint = vcat(degNbMaint, j-1)
            end
        end
    end  

    if (!isempty(degDates)) && (degDates[1] == 0.)
        degDates = degDates[2:end] # remove value of the degraddation at time 0
        degNbMaint = degNbMaint[2:end]
        degTypes = degTypes[2:end]
    end
    
    return DegradationData(
        DataFrame(DATE = maintenanceDates, TYPE = maintenanceTypes),
        DataFrame(DATE = degDates, NB_MAINTENANCES = degNbMaint, VALUE = fill(0.,length(degDates)), TYPE = degTypes)
        )
end

#####GENERATEUR DE DONNEES ALEATOIRES#####
function generateur(wienerARD1::WienerARD1, degradationdata::DegradationData)
    deg = degradationdata.degradations
    maint = degradationdata.maintenances

    μ = wienerARD1.underlyingDegradation[1]
    σ = wienerARD1.underlyingDegradation[2]
    type_to_ρ = Dict(wienerARD1.maintenances.TYPE .=> wienerARD1.maintenances.VALUE)

    NBM = filter(n -> n in deg.NB_MAINTENANCES, 0:nrow(maint)) #creates the vector ν
    ν = NBM .+ 1
    K = length(ν)

    #Y2 takes the successive simulated values of the degradation level over time
    #Y1 corresponds to Y(τ_j+) which is necessary to compute Y(τ_{j+1}+) = Y(τ_{j+1}-) - ρ[Y(τ_{j+1}-) - Y(τ_j+)] 
    Y1 = 0.
    Y2 = 0.

    new_deg = Float64[]#the vector which will contain the new values for degradationdataz.degradations.VALUE

    #The first increment is special 
    if ν[1] == 1
        Δt = deg[1, "DATE"]
        Y2 += rand(Normal(μ * Δt, σ * sqrt(Δt)))
        push!(new_deg, Y2)
    else
        Δτ = maint[1, "DATE"]
        Y2 += rand(Normal(μ * Δτ, σ * sqrt(Δτ)))
        Y1 = Y2 - type_to_ρ[maint[1, "TYPE"]] * (Y2 - Y1)

        if ν[1] != 2
            for j in 2:ν[1]-1
                Δτ = maint[j, "DATE"] - maint[j-1, "DATE"]
                Y2 = Y1 + rand(Normal(μ * Δτ, σ * sqrt(Δτ)))
                Y1 = Y2 - type_to_ρ[maint[j, "TYPE"]] * (Y2 - Y1)
            end
        end

        Δt = deg[1, "DATE"] - maint[deg[1, "NB_MAINTENANCES"], "DATE"]
        Y2 = Y1 + rand(Normal(μ * Δt, σ * sqrt(Δt)))
        push!(new_deg, Y2)
    end

    #For each ν[i], we generate y_{ν_i, j} for j≥2 and then generate y_{ν_{i+1}, 1} (if ν_{i+1} exists)
    for i in 1:K
        ν_ith_deg = filter(row -> row.NB_MAINTENANCES == ν[i] - 1, deg)#keep only the observations between \nu[i]-1 and \nu[i]

        if nrow(ν_ith_deg) != 1
            Δt = diff(ν_ith_deg.DATE)
            ν_ith_deg[2:end, "VALUE"] = new_deg[end] .+ cumsum(rand.(Normal.(μ .* Δt, σ .* sqrt.(Δt))))
        end

        if i != lastindex(ν)
            u_ν_i = maint[ν[i], "TYPE"]
            ρ_u_ν_i = type_to_ρ[u_ν_i]
            
            idx = findfirst(deg.NB_MAINTENANCES .== ν[i + 1] - 1)
            ν_iplus1th_first_deg = deg[[idx], :]#find the first observation following all the observations stored in ν_ith_deg

            Δt = maint[ν[i], "DATE"] - ν_ith_deg[end, "DATE"]
            Y2 = ν_ith_deg[end, "VALUE"] + rand(Normal(μ * Δt, σ * sqrt(Δt)))
            Y1 = Y2 - ρ_u_ν_i * (Y2 - Y1)

            if ν[i] + 1 != ν[i+1]
                for l in ν[i]+1:ν[i+1]-1
                    u_l = maint[l, "TYPE"]#N'y-aurait-il pas un l+1 au lieu de l ?
                    ρ_u_l = type_to_ρ[u_l]
                    Δτ_l = maint[l, "DATE"] - maint[l - 1, "DATE"]
                    Y2 = Y1 + rand(Normal(μ * Δτ_l, σ * sqrt(Δτ_l)))
                    Y1 = Y2 - ρ_u_l * (Y2 - Y1)
                end
            end

            Δt = ν_iplus1th_first_deg[1, "DATE"] - maint[ν[i], "DATE"]
            ν_iplus1th_first_deg[1, "VALUE"] = Y1 + rand(Normal(μ * Δt, σ * sqrt(Δt)))
            
            append!(new_deg, ν_ith_deg[2:end, "VALUE"])
            append!(new_deg, ν_iplus1th_first_deg.VALUE)
        else
            append!(new_deg, ν_ith_deg[2:end, "VALUE"])
        end
    end

    degradationdata.degradations.VALUE = new_deg

    return degradationdata
end