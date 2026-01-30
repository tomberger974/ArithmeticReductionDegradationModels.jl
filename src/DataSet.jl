"""
Data structure for observations on 1 degradation trajectory (real or simulated) and corresponding maintenances
"""
mutable struct DegradationData

    # Different caracteristics about the system from which the trajectory is issued
    # For real data that can allowed to select specific trajectories according to their caracteristics
    infos::Dict{String, Any}

    # Maintenances information
    # with 2 columns, called 
    #   DATE (Float64), containing maintenance times,
    #   TYPE (String), specifying maintenance type (not used now ...)
    maintenances::DataFrame
    
    #All degradation information, at least 3 columns:
    #   DATE (Float64): time at which degradation is observed,
    #   NB_MAINTENANCES (Int64): cummulative number of maintenance that have been done before degradation observation,
    #   VALUE (Float64): observed degradation value,
    #other columns can be added, in particular 
    #   TYPE (String): an indicator to specify different type of degradation observations
    degradations::DataFrame

    
    """
    DegradationData(
        maintenances::DataFrame, 
        degradations::DataFrame,
        infos::Dict{String, Any}=Dict{String, Any}())

    The constructor for degradation (and corresponding maintenancs) data set object where
        `maintenances` : a DataFrame with 2 columns, named `DATE` and `TYPE`, containing the successive maintenance times and their types
        `degradation` : a DataFrame with at least 3 columns, named `DATE`, `NB_MAINTENANCES`, `VALUE`, and potentially an additionnal one named `TYPE`, containing the successive degradation observation times, the number of maintenances already done before the degradation observation, the corresponding observed values of the degradation, and potentially type indicators of about the observation
        `infos` : some potential additionnal informations about the data set
    """
    function DegradationData(
        maintenances::DataFrame, 
        degradations::DataFrame,
        infos::Dict{String, T} where T<: Any=Dict{String, Any}())
        res = verifDegradationData(maintenances, degradations)
        new(convert(Dict{String, Any}, infos), res[1], res[2])
    end
end

#Verify that Dataframes maintenances and degradations respect the imposed format and are coherent
function verifDegradationData(maintenances::DataFrame, degradations::DataFrame)
    if sum([i ∉ names(degradations) for i ∈ ["DATE", "NB_MAINTENANCES", "VALUE"]]) > 0
        error("degradations DataFrame must contain columns DATE, NB_MAINTENANCES and VALUE")
    end
    if typeof(degradations.DATE) != Vector{Float64}
        degradations.DATE = convert(Vector{Float64}, degradations.DATE)
    end
    if typeof(degradations.VALUE) != Vector{Float64}
        degradations.VALUE = convert(Vector{Float64}, degradations.VALUE)
    end
    if typeof(degradations.NB_MAINTENANCES) != Vector{Int64}
        degradations.NB_MAINTENANCES = convert(Vector{Int64}, degradations.NB_MAINTENANCES)
    end
    if "DATE" ∉ names(maintenances)
        error("maintenances DataFrame must contain columns DATE")
    end
    if typeof(maintenances.DATE) != Vector{Float64}
        maintenances.DATE = convert(Vector{Float64}, maintenances.DATE)
    end
    if "TYPE" ∉ names(maintenances)
        maintenances[!, :TYPE] .= ""
    end
    if typeof(maintenances.TYPE) != Vector{String}
        error("maintenance types must be repesented with String")
    end

    #if necessary sort Dataframes values according to time occurence and previous maintenance observations
    if !issorted(maintenances.DATE)
        isort = sortperm(maintenances.DATE)
        maintenances = maintenances[isort,:]
    end
    if (!isempty(degradations.NB_MAINTENANCES)) && (maximum(degradations.NB_MAINTENANCES) > nrow(maintenances))
        error("The row NB_MAINTENANCES of DataFrame degradation is not coherent with the number of maintenance action in DataFrame maintenances")
    end
    if (!issorted(degradations.DATE)) .|| (!issorted(degradations.NB_MAINTENANCES))
        isort = sortperm(degradations.NB_MAINTENANCES)
        degradations = degradations[isort,:]
        inb = 0
        for i in 0:maximum(degradations.NB_MAINTENANCES)
            deg = degradations[degradations.NB_MAINTENANCES .== i,"DATE"]
            if (i>0) && (sum(deg .< maintenances.DATE[i]) > 0)
                error("Incoherence(s) between the column NB_MAINTENANCES of degradation and DATE of maintenances and degradations")
            end
            isort[ (inb+1): (inb+length(deg)) ] =  inb .+ sortperm(deg)
            inb += length(deg)
        end
        degradations = degradations[isort,:]
        if (!issorted(degradations.DATE)) .|| (!issorted(degradations.NB_MAINTENANCES))
            error("Incoherence(s) between column DATE and NB_MAINTENANCES degradation DataFrame")
        end
    end 
    return (maintenances, degradations)
end

function DegradationData(maintenanceDates::Vector{Float64}, Δt::Float64, t_max::Float64=maintenanceDates[end], maintenanceTypes::Vector{String}=fill("",length(maintenanceDates)); imaint_from::Int64=0, t_from::Float64=0.)

    #to avoid errors
    if Δt <= 0.
        error("Time increment Δt must be positiv")
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

"""
    maintenances(data::DegradationData)

    Return a copy of maintenances Dataframe of the `data` set
"""
function maintenances(data::DegradationData)
    return deepcopy(data.maintenances)
end

"""
    degradations(data::DegradationData)

    Return a copy of degradations DataFrame of the `data` set
"""
function degradations(data::DegradationData)
    return deepcopy(data.degradations)
end

"""
    maintenances!(data::DegradationData, maintenances::DataFrame)

    Change maintenances of the `data` set. For format constraints on `maintenances` see function `DegradationData`
"""
function maintenances!(data::DegradationData, maintenances::DataFrame)
    res = verifDegradationData(maintenances, degradations(data))
    data.maintenances = res[1]
    return data
end

"""
    degradations!(data::DegradationData, degradations::DataFrame)

    Change degradations of the `data` set. For format constraints on `degradations` see function `DegradationData`
"""
function degradations!(data::DegradationData, degradations::DataFrame)
    res = verifDegradationData(maintenances(data), degradations)
    data.degradations = res[2]
    return data
end

"""
    degradationsANDmaintenances!(data::DegradationData, maintenances::DataFrame, degradations::DataFrame)

    Change degradations and maintenances of the `data` set. For format constraints on `degradations` and maintenances see function `DegradationData`
"""
function degradationsANDmaintenances!(data::DegradationData, maintenances::DataFrame, degradations::DataFrame)
    res = verifDegradationData(maintenances, degradations)
    data.maintenances = res[1]
    data.degradations = res[2]
    return data
end

"""
    infos(data::DegradationData)

    Return  a copy of `infos` Dict of `data` set
"""
function infos(data::DegradationData)
    return deepcopy(data.infos)
end

"""
    infos!(data::DegradationData, infos::Dict{String, Any})

    Modify optionnal informations in `data`
"""
function infos!(data::DegradationData, infos::Dict{String, T} where T<: Any=Dict{String, Any}())
    data.infos = infos
    return data
end