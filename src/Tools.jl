function min_max(x::Vector{Float64}, y::Vector{Float64})
    if !isempty(y)
        return [min(x[1], minimum(y)), max(x[2], maximum(y))]
    else
        return x
    end
end

function min_max(x::Vector{Float64}, y::Vector{Any})
    if !isempty(y)
        return [min(x[1], minimum(y)), max(x[2], maximum(y))]
    else
        return x
    end
end

# Internal function which provides the index of a given maintenance type 
function which_maint_type(type::String, maintTypes::Vector{String})
    if isempty(maintTypes)
        return 1
    else
        for (i,mtype) in enumerate(maintTypes)
            if mtype == type
                return i
            end
        end
    end
    @warn string("The maintenance type ", type, " apear in the data and has not been defined, consequantly it is assimilated to type ", maintTypes[1])
    return 1
end

"""
    CreateData(maintenancesDates::Vector{Float64}, Δt::Float64, t_max::Float64=maintenancesDates[end], maintenancesTypes::Vector{Float64}=fill("",length(maintenancesDates)))

    Create a degradation data set `DegradationData` (observed values of the degradation are undef)
        - maintenances at times `maintenancesDates` with types that can be specifed by `maintenancesTypes`
        - periodic degradation observations every `Δt` units of time 
        - from 0 or if specified the max between the `imaint_from`-th maintenance time and time `t_from`
        - up to the min between the last maintenance time and `t_max` if it is specified
    Degradation observations have a TYPE specifying if there are done "After", "Before" or "Between" maintenances

"""
function CreateData(maintenanceDates::Vector{Float64}, Δt::Float64, t_max::Float64=maintenanceDates[end], maintenanceTypes::Vector{String}=fill("",length(maintenanceDates)); imaint_from::Int64=0, t_from::Float64=0.)

    if Δt <= 0.
        error("Time increment Δt must be positiv")
    end
    if length(maintenanceDates) != length(maintenanceTypes)
        error("Vectors of maintenancesDates and maintenancesTypes must have the same length.")
    end

    degDates = Vector{Float64}()
    degNbMaint = Vector{Int64}()
    degTypes = Vector{String}()

    if imaint_from == 0
        maint = vcat(0,maintenanceDates[1:end])
    else
        maint = maintenanceDates[imaint_from:end]
    end

    for (jj,τ) in enumerate(maint)
        if τ > t_max 
            break
        end
        j = jj + imaint_from
        if (τ >= t_from)
            degDates = vcat(degDates, τ)
            degNbMaint = vcat(degNbMaint, j-1)
            degTypes = vcat(degTypes, "After")
        end
        t = max(t_from, τ + Δt)
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