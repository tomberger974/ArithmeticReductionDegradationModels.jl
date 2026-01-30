abstract type MaintenancesPlot end
abstract type DegradationsPlot end
abstract type QuantilesPlot end

import CairoMakie.plot!

"""
    plot!(ax::Axis, data::DegradationData, ::Type{MaintenancesPlot}, type::Union{Missing, String}=missing; color::Any = :black, linestyle::Union{Symbol,Char}= :dash, optargs...)

    On an Axis object `ax` define on a Figure object, add the representation of maintenance times, that is to say,
    vertical lines at the successive maintenance times of a `data` DegradationData data set (column `DATE` of the DataFrame `data.maintenances`)
    Optionnaly it is possible to plot lines only at maintenance times of certain `type`, i.e. the ones for which the colum `TYPE` of the DataFrame `data.maintenances` corresponds to the value of `type`.
    Optional arguments are passed to `vlines`, for exemple color or linestyle arguments

"""
function plot!(ax::Axis, data::DegradationData, ::Type{MaintenancesPlot}, type::Union{Missing, String}=missing; color::Any = :black, linestyle::Union{Symbol,Char}= :dash, optargs...)
    if (linestyle != :none) 
        for i in 1:length(data.maintenances.DATE)
            if ismissing(type) || (data.maintenances.TYPE[i] == type)
                vlines!(ax, data.maintenances.DATE[i]; color=color, linestyle=linestyle, optargs...)
            end
        end
    end
end

"""
    plot_deg!(ax::Axis, data::DegradationData, ::Type{DegradationPlot}, types::Vector{String}=Vector{String}(); 
            marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
            linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash,  
            plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    On an Axis object `ax` define on a Figure object, add the representation of the degradation data, that is to say,
    Plot the successive degradation values versus time of a `data` DegradationData data set (column `DATE` versus column `VALUE` of the DataFrame `data.degradations`).
    By default degradation observations are connected by a line when they take place between the same 2 successive maintenances, use optional argument `linestyle = :none` to remove or modify the style of those lines.
    By default degradation observations are marked with a special diamond maker which is complet for observations done bewteen maintence actions and left-half or right-half for observations respectively done just before and just after maintenance actions. Marker can be removed with optional argument `marker = :none`, or replace by any other marker types, like for exemple `marker='x'` for using crosses.
    Line and marker colors can be specified with respectivelly optional arguments `linecolor` and `markercolor`.
    Optional argument `plot_cens = true`, can be used in the case where degradation `data` have a supplementary column named `CENS`, in order to additionnaly plot a special up arrow marker on all degradation observations having value 1 in this column.
    Optional argument `plot_pb_nb_maint = true`, can be used in the case where degradation `data` have a supplementary column named `PB_NB_MAINTENANCES`, in order to plot a complet diamond on the observation done at maintenance times (instead of half-left or half-right diamond depending on the fact that obsrevation takes place befor or after maintenance), but this is only for degradation observation times forwhich there is the value 1 in this column.
"""
function plot!(ax::Axis, data::DegradationData, ::Type{DegradationsPlot}, types::Vector{String}=Vector{String}(); 
            marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
            linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash,  
            plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    if plot_cens && ("CENS" ∉ names(data.degradations))
        error("Optional argument 'plot_cens = true' needs to specify a column 'CENS' in the degradation data.")
    end
    if plot_pb_nb_maint && ("PB_NB_MAINTENANCES" ∉ names(data.degradations))
        error("Optional argument 'plot_pb_nb_maint = true' needs to specify a column 'PB_NB_MAINTENANCES' in the degradation data.")
    end
    
    #the special up arrow marker to optionnaly represent censored degradation data
    arrow_path = BezierPath([
        MoveTo(Point(0, 0)),
        LineTo(Point(0.05, 0)),
        LineTo(Point(0.05, 1.25)),
        LineTo(Point(0.6, .75)),
        LineTo(Point(0, 1.5)),
        LineTo(Point(-0.6, 0.75)),
        LineTo(Point(-0.05, 1.25)),
        LineTo(Point(-0.05, 0)),
        ClosePath()
    ])

    #special marker of half-right and half-left diamond for degradation observations just befor or just after maintenance
    rtriangle_path = BezierPath([
        MoveTo(Point(0, 0)),
        LineTo(Point(0, 0.5)),
        LineTo(Point(0.5, 0)),
        LineTo(Point(0, -.5)),
        ClosePath()
    ])
    ltriangle_path = BezierPath([
        MoveTo(Point(0, 0)),
        LineTo(Point(0, 0.5)),
        LineTo(Point(-0.5, 0)),
        LineTo(Point(0, -.5)),
        ClosePath()
    ])
    #and the corresponding complet diamond marker for degradation observations between maintenances
    diamond_path = BezierPath([
        MoveTo(Point(0, 0.5)),
        LineTo(Point(0.5, 0)),
        LineTo(Point(0, -0.5)),
        LineTo(Point(-0.5, 0)),
        ClosePath()
    ])

    minmaxy = [Inf64, -Inf64]#Minimal and maximal degradation values, in order to potentialy choose the bound on the ordinate axis
    prev_maint = 0
    if !isempty(data.degradations)
        for i in 0:maximum(data.degradations.NB_MAINTENANCES)
            if (length(data.maintenances.DATE) > i)
                suiv_maint = data.maintenances.DATE[i+1]#next maintenance time
            else
                suiv_maint = Inf64#we are currently after the last maintenance time
            end 
            if (length(types) == 0)#do we represent all degradation observations or only those of types `types``
                deg = data.degradations[
                    (data.degradations.NB_MAINTENANCES .== i)  , #we consider only degradation data done after the ith maintenance and before the (i+1)th maintenance
                    :]
            else
                deg = data.degradations[
                    (data.degradations.NB_MAINTENANCES .== i) .&& #we consider only degradation data done after the ith maintenance and before the (i+1)th maintenance
                    [x  ∈ types for x in data.degradations.TYPE] , 
                    :]
            end

            minmaxy = min_max(minmaxy, deg.VALUE) #minmax function defined in Tools.jl
            if !isempty(deg)
                if linestyle != :none
                    lines!(ax, deg.DATE, deg.VALUE; color=linecolor, linestyle = linestyle)
                end
                if marker == :deg
                    degm = deg[(deg.DATE .> prev_maint) .&& (deg.DATE .< suiv_maint),:]#degradation observation done between maintenances
                    if plot_pb_nb_maint
                        degpb = deg[((deg.DATE .== prev_maint) .|| (deg.DATE .== suiv_maint)) .&& (deg.PB_NB_MAINTENANCES .== 1), :]#degradation observations done in maintenance the ith or (i+1)th maintenance time and for which the PB_NB_MAINTENANCES column has a value equal to 1
                        degl = deg[(deg.DATE .== prev_maint) .&& (deg.PB_NB_MAINTENANCES .!= 1), :]
                        degr = deg[(deg.DATE .== suiv_maint) .&& (deg.PB_NB_MAINTENANCES .!= 1), :]
                    else
                        degpb = []
                        degl = deg[deg.DATE .== prev_maint, :]# degradation observations done just after the ith maintenance
                        degr = deg[deg.DATE .== suiv_maint, :]#degradation observationd done just before the (i+1)th maintenance
                    end
                    if !isempty(degl)
                        Makie.scatter!(ax, degl.DATE, degl.VALUE; color=markercolor, marker=rtriangle_path, markersize = markersize)
                    end
                    if !isempty(degpb)
                        Makie.scatter!(ax, degpb.DATE, degpb.VALUE; color=markercolor, marker=diamond_path, markersize = markersize)
                    end
                    if !isempty(degm)
                        Makie.scatter!(ax, degm.DATE, degm.VALUE; color=markercolor, marker=diamond_path, markersize = markersize)
                    end
                    if !isempty(degr)
                        Makie.scatter!(ax, degr.DATE, degr.VALUE; color=markercolor, marker=ltriangle_path, markersize = markersize)
                    end 
                else
                    if marker != :none
                        Makie.scatter!(ax, deg.DATE, deg.VALUE; color=markercolor, marker=marker, markersize = markersize)
                    end
                end
                if plot_cens
                    degc = deg[deg.CENS .== 1, :]
                    if !isempty(degc)
                        Makie.scatter!(ax, degc.DATE, degc.VALUE; color=markercolor, marker=arrow_path, markersize = markersize)
                    end
                end
            end
            prev_maint = suiv_maint
        end
    end
    maxx = 0
    if !isempty(data.degradations)
        maxx = max(maxx, data.degradations.DATE[end])
    end
    if !isempty(data.maintenances)
        maxx = max(maxx, data.maintenances.DATE[end])
    end
    return ( [0, maxx], minmaxy)
end

"""
    plot!(ax::Axis, model::WienerARD∞, data::DegradationData, ::Type{QuantilesPlot}, probs::Vector{Float64}=[0.25, 0.5, 0.75]; 
    color::Any=:red, colors:: Vector{T} where T <: Any = fill(color, length(probs)), linestyles::Union{Vector{Symbol},Vector{Char}}= fill(:dot, length(probs)),    
    froms::Vector{Int64}=[0], tmaxs::Vector{Float64}= [max(data.maintenances.DATE[end], data.degradations.DATE[end])], step::Int64 = 500,
    withband::Bool=true, colorband::Any=(color, 0.3) )

    On an Axis object `ax` define on a Figure object, add the representation quantiles of the degradation at successive times according the the WienerARD∞ model `model`,
     - `probs` is a vector of probabilies for which quantiles are represented
     - `step` indicates the step on abscissa line for quantile representation
     - quantiles are computed given all the degradation observations up to the `froms`-th one (for each element of `froms`)
     - quantiles are represented uf to time `t_maxs` (for each element of `t_maxs`)
     CARREFUL `t_maxs` and `from`s must have the same length
     - if `withband=true` a band is also represented between quantils of level mnimum(`probs`) and maximum(probs)
     - `color` allows to speficied a unique color used for lines and band representation, otherwise the band color can be specified with `colorband` and the lines colors with the vector `colors`
     - the lines styles can be specified with Vecor optionnal argument `linestyles`
"""
function plot!(ax::Axis, model::WienerARD∞, data::DegradationData, ::Type{QuantilesPlot}, probs::Vector{Float64}=[0.25, 0.5, 0.75]; 
    color::Any=:red, colors:: Vector{T} where T <: Any = fill(color, length(probs)), linestyles::Union{Vector{Symbol},Vector{Char}}= fill(:dot, length(probs)),    
    froms::Vector{Int64}=[0], tmaxs::Vector{Float64}= [max(data.maintenances.DATE[end], data.degradations.DATE[end])], step::Int64 = 500,
    withband::Bool=true, bandcolor::Any=(color, 0.3) )

    if length(probs) != length(colors)
        error("Vector arguments linecolors and probs must have the same length")
    end
    if length(probs) != length(linestyles)
        error("Vector arguments linecolors and probs must have the same length")
    end
    if length(froms) != length(tmaxs)
        error("Arguments from and tmax must have the same length")
    end
    if minimum(froms) < 0
        error("Element(s) of argument from must be positiv")
    end
    if maximum(froms)  > nrow(data.degradations)
        error("Element(s) of rgument from can nnot be greater than degradation observation number")
    end

    minmaxy = [0, -Inf64]
    maxx = -Inf64

    for i in 1:length(froms)
        from = froms[i]
        tmax = tmaxs[i]
        println((from, tmax))

        if from == 0 #not conditional, we only know that degradation is equals to 0 at time 0
            qs = CreateData(data.maintenances.DATE, convert(Float64, tmax/step), tmax, data.maintenances.TYPE)# function CreateData in Tools.jl
            prepend!(qs.degradations, DataFrame(DATE=[0.], NB_MAINTENANCES= [0], VALUE = [0.], TYPE = "Beween"))
            qss = quantile(model, qs, probs)
            for nbmaint in 0:nrow(data.maintenances)
                qs_nbmaint = qss[qss.NB_MAINTENANCES .== nbmaint, :]
                if withband
                    Makie.band!(ax, qs_nbmaint.DATE, qs_nbmaint[!,string("Q_", minimum(probs))], qs_nbmaint[!,string("Q_", maximum(probs))], color=bandcolor)
                end
                for (j, prob) in enumerate(probs)
                    lines!(ax, qs_nbmaint.DATE, qs_nbmaint[!,string("Q_", prob)], color = colors[j], linestyle=linestyles[j])
                end
            end
        else
            qs = CreateData(data.maintenances.DATE, convert(Float64, tmax/step), tmax, data.maintenances.TYPE, imaint_from = data.degradations.NB_MAINTENANCES[from], t_from = data.degradations.DATE[from]) # function CreateData in Tools.jl
            if !isempty(qs.degradations.VALUE)
                qs.degradations.VALUE[1] = data.degradations.VALUE[from]
                qss = quantile(model, qs, probs, from=1)
                if (qs.degradations.TYPE == "After") && (nrow(qss) > 1)
                    qss = qss[1:(nrow(qss)-1), :] #if plot ends in a maintenance time, we do not plot quantiles after maintenance
                end
                for nbmaint in unique(qss.NB_MAINTENANCES)
                    qs_nbmaint = qss[qss.NB_MAINTENANCES .== nbmaint, :]
                    for (j, prob) in enumerate(probs)
                        lines!(ax, qs_nbmaint.DATE, qs_nbmaint[!,string("Q_", prob)], color = colors[j], linestyle=linestyles[j])
                    end
                    if withband
                        Makie.band!(ax, qs_nbmaint.DATE, qs_nbmaint[!,string("Q_", minimum(probs))], qs_nbmaint[!,string("Q_", maximum(probs))], color = bandcolor)
                    end
                end
            end
        end
        if !isempty(qs.degradations.VALUE)
            minmaxy[1] = min(minmaxy[1], minimum(qss[!,string("Q_", minimum(probs))])) 
            minmaxy[2] = min(minmaxy[2], maximum(qss[!,string("Q_", maximum(probs))]))
            maxx = max(maxx,  qss.DATE[nrow(qss)])
        end
    end

    return ( [0, maxx], minmaxy)
end

"""
    plot!(ax::Axis, data::DegradationData, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    On an Axis object `ax` define on a Figure object, add the representation of the degradation and maintenance data, that is to say,
    Plot the successive degradation values versus time of a `data` DegradationData data set (column `DATE` versus column `VALUE` of the DataFrame `data.degradations`).
    and plot also vertical lines at the successive maintenance times of a `data` DegradationData data set (column `DATE` of the DataFrame `data.maintenances`)

    Optional arguments for degradation observations :
     - By default degradation observations are connected by a line when they take place between the same 2 successive maintenances, use optional argument `linestyle = :none` to remove or modify the style of those lines.
     - By default degradation observations are marked with a special diamond maker which is complet for observations done bewteen maintence actions and left-half or right-half for observations respectively done just before and just after maintenance actions. Marker can be removed with optional argument `marker = :none`, or replace by any other marker types, like for exemple `marker='x'` for using crosses.
     - Line and marker colors representing degradation observations can be specified with respectivelly optional arguments `linecolor` and `markercolor`.
     - Optional argument `plot_cens = true`, can be used in the case where degradation `data` have a supplementary column named `CENS`, in order to additionnaly plot a special up arrow marker on all degradation observations having value 1 in this column.
     - Optional argument `plot_pb_nb_maint = true`, can be used in the case where degradation `data` have a supplementary column named `PB_NB_MAINTENANCES`, in order to plot a complet diamond on the observation done at maintenance times (instead of half-left or half-right diamond depending on the fact that obsrevation takes place befor or after maintenance), but this is only for degradation observation times forwhich there is the value 1 in this column.
    
     The color and linestyle of vertical lines representing maintenance actions can be specified by optional arguments:
     - `maintenancescolors` a DataFrame with 2 columns named `TYPE` and `COLOR` whose lines must contains all the maintenance types and the corresponding choosen color,
     - `maintenanceslinestyles` a DataFrame with 2 columns named `TYPE` and `LINESTYLE` whose lines must contains all the maintenance types and the corresponding choosen linestyle,
     - `maintenancecolor`, if `maintenancescolors` argument is not specified, this is a unique color argument for all types oof maintenances,
     - `maintenancelinestyle`, if `maintenanceslinestyles` argument is not specified, this is a unique linestyle argument for all types oof maintenances,
"""
function plot!(ax::Axis, data::DegradationData, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    if isempty(maintenancescolors)
        maintenancescolors = DataFrame(TYPE = unique(data.maintenances.TYPE), COLOR = fill(maintenancecolor, length(unique(data.maintenances.TYPE))))
    else
        if sum([i ∉ names(maintenancescolors) for i ∈ ["TYPE", "COLOR"]]) > 0
            error("maintenancescolors DataFrame must contain columns TYPE and COLOR")
        end
        if !isempty(data.maintenances) && (sum([i ∉ unique(maintenancescolors.TYPE) for i ∈ unique(data.maintenances.TYPE)]) > 0)
            error("Column TYPE of maintenancescolors DataFrame must contain all types of maintenances")
        end
    end
    if isempty(maintenanceslinestyles)
        maintenanceslinestyles = DataFrame(TYPE = unique(data.maintenances.TYPE), LINESTYLE = fill(maintenancelinestyle, length(unique(data.maintenances.TYPE))))
    else
        if sum([i ∉ names(maintenanceslinestyles) for i ∈ ["TYPE", "LINESTYLE"]]) > 0
            error("maintenanceslinestyles DataFrame must contain columns TYPE and LINESTYLE")
        end
        if !isempty(data.maintenances) && (sum([i ∉ unique(maintenanceslinestyles.TYPE) for i ∈ unique(data.maintenances.TYPE)]) > 0)
            error("Column TYPE of maintenanceslinestyles DataFrame must contain all types of maintenances")
        end
    end

    for type in unique(data.maintenances.TYPE)
        plot!(ax, data, MaintenancesPlot, type ; color = maintenancescolors.COLOR[maintenancescolors.TYPE .== type][1], linestyle = maintenanceslinestyles.LINESTYLE[maintenanceslinestyles.TYPE .== type][1])
    end
    
    bound = plot!(ax, data, DegradationsPlot, degtypes; marker=marker, markersize=markersize, markercolor=markercolor, linecolor=linecolor, linestyle=linestyle, plot_cens=plot_cens, plot_pb_nb_maint=plot_pb_nb_maint)
    if marker == :deg
        plot!(ax, data, DegradationsPlot, degtypes; marker=:circle, markersize=round(Int64, markersize/4), markercolor=:white, linecolor=linecolor, linestyle=linestyle, plot_cens=plot_cens, plot_pb_nb_maint=plot_pb_nb_maint)
    end
    return bound
end

"""
    plot!(f::Figure, data::DegradationData, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    xlims::Union{Missing, Vector{Float64}}=missing, ylims::Union{Missing, Vector{Float64}}=missing,
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    On a Figure object `f`, define an Axis object in order to represent the degradation and maintenance data, that is to say,
    Plot the successive degradation values versus time of a `data` DegradationData data set (column `DATE` versus column `VALUE` of the DataFrame `data.degradations`).
    and plot also vertical lines at the successive maintenance times of a `data` DegradationData data set (column `DATE` of the DataFrame `data.maintenances`)

    Optional arguments for degradation observations :
     - By default degradation observations are connected by a line when they take place between the same 2 successive maintenances, use optional argument `linestyle = :none` to remove or modify the style of those lines.
     - By default degradation observations are marked with a special diamond maker which is complet for observations done bewteen maintence actions and left-half or right-half for observations respectively done just before and just after maintenance actions. Marker can be removed with optional argument `marker = :none`, or replace by any other marker types, like for exemple `marker='x'` for using crosses.
     - Line and marker colors representing degradation observations can be specified with respectivelly optional arguments `linecolor` and `markercolor`.
     - Optional argument `plot_cens = true`, can be used in the case where degradation `data` have a supplementary column named `CENS`, in order to additionnaly plot a special up arrow marker on all degradation observations having value 1 in this column.
     - Optional argument `plot_pb_nb_maint = true`, can be used in the case where degradation `data` have a supplementary column named `PB_NB_MAINTENANCES`, in order to plot a complet diamond on the observation done at maintenance times (instead of half-left or half-right diamond depending on the fact that obsrevation takes place befor or after maintenance), but this is only for degradation observation times forwhich there is the value 1 in this column.
    
     The color and linestyle of vertical lines representing maintenance actions can be specified by optional arguments:
     - `maintenancescolors` a DataFrame with 2 columns named `TYPE` and `COLOR` whose lines must contains all the maintenance types and the corresponding choosen color,
     - `maintenanceslinestyles` a DataFrame with 2 columns named `TYPE` and `LINESTYLE` whose lines must contains all the maintenance types and the corresponding choosen linestyle,
     - `maintenancecolor`, if `maintenancescolors` argument is not specified, this is a unique color argument for all types oof maintenances,
     - `maintenancelinestyle`, if `maintenanceslinestyles` argument is not specified, this is a unique linestyle argument for all types oof maintenances,

"""
function plot!(f::Figure, data::DegradationData, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    xlims::Union{Missing, Vector{Float64}}=missing, ylims::Union{Missing, Vector{Float64}}=missing,
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    ax = Axis(f[1,1])
    bounds = plot!(ax, data, degtypes; marker=marker, markersize=markersize, markercolor=markercolor, linecolor=linecolor, linestyle=linestyle, maintenancecolor=maintenancecolor, maintenancelinestyle=maintenancelinestyle, maintenancescolors=maintenancescolors, maintenanceslinestyles=maintenanceslinestyles, plot_cens=plot_cens, plot_pb_nb_maint=plot_pb_nb_maint)
    if ismissing(xlims)
        xlims = bounds[1]
    end
    if ismissing(ylims)
        ylims = bounds[2]
    end
    Makie.xlims!(ax, 0, xlims[2]*1.05)
    Makie.ylims!(ax, min(0, ylims[1]), ylims[2]*1.05)
    
    return (ax, bounds)
end

"""
    plot!(f::Figure, datas::Vector{DegradationData}, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    xlims::Union{Missing, Vector{Float64}}=missing, ylims::Union{Missing, Vector{Float64}}=missing,
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false, nbcols::Int64=1)

    On a Figure object `f`, define an Axis object for each DegradationData element of the vector `data` in order to represent the corresponding degradation and maintenance data, that is to say,
    Plot the successive degradation values versus time of a `data` DegradationData data set (column `DATE` versus column `VALUE` of the DataFrame `data.degradations`).
    and plot also vertical lines at the successive maintenance times of a `data` DegradationData data set (column `DATE` of the DataFrame `data.maintenances`)

    Optional argument nbcols specifies how the different representation are organized on the Figure by specying the number of plot per column
    If an element of `datas` has in its `infos` an element with key `name` its value is use as the title of the corresponding plot, otherwise it has no title

    Optional arguments for degradation observations :
     - By default degradation observations are connected by a line when they take place between the same 2 successive maintenances, use optional argument `linestyle = :none` to remove or modify the style of those lines.
     - By default degradation observations are marked with a special diamond maker which is complet for observations done bewteen maintence actions and left-half or right-half for observations respectively done just before and just after maintenance actions. Marker can be removed with optional argument `marker = :none`, or replace by any other marker types, like for exemple `marker='x'` for using crosses.
     - Line and marker colors representing degradation observations can be specified with respectivelly optional arguments `linecolor` and `markercolor`.
     - Optional argument `plot_cens = true`, can be used in the case where degradation `data` have a supplementary column named `CENS`, in order to additionnaly plot a special up arrow marker on all degradation observations having value 1 in this column.
     - Optional argument `plot_pb_nb_maint = true`, can be used in the case where degradation `data` have a supplementary column named `PB_NB_MAINTENANCES`, in order to plot a complet diamond on the observation done at maintenance times (instead of half-left or half-right diamond depending on the fact that obsrevation takes place befor or after maintenance), but this is only for degradation observation times forwhich there is the value 1 in this column.
    
     The color and linestyle of vertical lines representing maintenance actions can be specified by optional arguments:
     - `maintenancescolors` a DataFrame with 2 columns named `TYPE` and `COLOR` whose lines must contains all the maintenance types and the corresponding choosen color,
     - `maintenanceslinestyles` a DataFrame with 2 columns named `TYPE` and `LINESTYLE` whose lines must contains all the maintenance types and the corresponding choosen linestyle,
     - `maintenancecolor`, if `maintenancescolors` argument is not specified, this is a unique color argument for all types oof maintenances,
     - `maintenancelinestyle`, if `maintenanceslinestyles` argument is not specified, this is a unique linestyle argument for all types oof maintenances,

"""
function plot!(f::Figure, datas::Vector{DegradationData}, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    xlims::Union{Missing, Vector{Float64}}=missing, ylims::Union{Missing, Vector{Float64}}=missing,
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false, nbcols::Int64=1)

    axs = Vector{Axis}(undef,length(datas))
    for i in 1:length(datas)
        ligne = floor(Int, (i-1)/nbcols)
        colonne = i - ligne * nbcols
        if haskey(datas[i].infos, "name")
            axs[i] = Axis(f[ligne,colonne], title = string(datas[i].infos["name"]))
        else
            axs[i] = Axis(f[ligne,colonne])
        end
    end
    bounds = plot!(axs, datas, degtypes; marker=marker, markersize=markersize, markercolor=markercolor, linecolor=linecolor, linestyle=linestyle, maintenancecolor=maintenancecolor, maintenancelinestyle=maintenancelinestyle, maintenancescolors=maintenancescolors, maintenanceslinestyles=maintenanceslinestyles, xlims=xlims, ylims=ylims, plot_cens=plot_cens, plot_pb_nb_maint=plot_pb_nb_maint)
    
    return (axs, bounds)
end

"""
    plot!(axs::Vector{Axis}, datas::Vector{DegradationData}, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    xlims::Union{Missing, Vector{Float64}}=missing, ylims::Union{Missing, Vector{Float64}}=missing,
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    On a Vector of Axis object add the representations corresponding to the different DegradationData element of the vector `data` in order to represent the corresponding degradation and maintenance data, that is to say,
    Plot the successive degradation values versus time of a `data` DegradationData data set (column `DATE` versus column `VALUE` of the DataFrame `data.degradations`).
    and plot also vertical lines at the successive maintenance times of a `data` DegradationData data set (column `DATE` of the DataFrame `data.maintenances`).
    Be carreful, `axs` and `datas` must have the same length.

    Optional arguments for degradation observations :
     - By default degradation observations are connected by a line when they take place between the same 2 successive maintenances, use optional argument `linestyle = :none` to remove or modify the style of those lines.
     - By default degradation observations are marked with a special diamond maker which is complet for observations done bewteen maintence actions and left-half or right-half for observations respectively done just before and just after maintenance actions. Marker can be removed with optional argument `marker = :none`, or replace by any other marker types, like for exemple `marker='x'` for using crosses.
     - Line and marker colors representing degradation observations can be specified with respectivelly optional arguments `linecolor` and `markercolor`.
     - Optional argument `plot_cens = true`, can be used in the case where degradation `data` have a supplementary column named `CENS`, in order to additionnaly plot a special up arrow marker on all degradation observations having value 1 in this column.
     - Optional argument `plot_pb_nb_maint = true`, can be used in the case where degradation `data` have a supplementary column named `PB_NB_MAINTENANCES`, in order to plot a complet diamond on the observation done at maintenance times (instead of half-left or half-right diamond depending on the fact that obsrevation takes place befor or after maintenance), but this is only for degradation observation times forwhich there is the value 1 in this column.
    
     The color and linestyle of vertical lines representing maintenance actions can be specified by optional arguments:
     - `maintenancescolors` a DataFrame with 2 columns named `TYPE` and `COLOR` whose lines must contains all the maintenance types and the corresponding choosen color,
     - `maintenanceslinestyles` a DataFrame with 2 columns named `TYPE` and `LINESTYLE` whose lines must contains all the maintenance types and the corresponding choosen linestyle,
     - `maintenancecolor`, if `maintenancescolors` argument is not specified, this is a unique color argument for all types oof maintenances,
     - `maintenancelinestyle`, if `maintenanceslinestyles` argument is not specified, this is a unique linestyle argument for all types oof maintenances,

"""
function plot!(axs::Vector{Axis}, datas::Vector{DegradationData}, degtypes::Vector{String}=Vector{String}(); 
    marker::Union{Symbol,Char}=:deg, markersize::Int64=20, markercolor::Any=:black,
    linecolor::Any=:purple, linestyle::Union{Symbol,Char}=:dash, 
    maintenancecolor::Any=:black, maintenancelinestyle::Union{Symbol,Char}=:dot,
    maintenancescolors::DataFrame=DataFrame(), maintenanceslinestyles::DataFrame=DataFrame(),
    xlims::Union{Missing, Vector{Float64}}=missing, ylims::Union{Missing, Vector{Float64}}=missing,
    plot_cens::Bool=false, plot_pb_nb_maint::Bool=false)

    if length(axs) != length(datas)
        error("Axis vector axs and DagradationData vector datas must have the same length")
    end
    minmaxy = [Inf64, -Inf64]#Minimal and maximal degradation values, in order to potentialy choose the bound on the ordinate axis
    maxx = 0
    for i in 1:length(datas)
        bounds = plot!(axs[i], datas[i], degtypes; marker=marker, markersize=markersize, markercolor=markercolor, linecolor=linecolor, linestyle=linestyle, maintenancecolor=maintenancecolor, maintenancelinestyle=maintenancelinestyle, maintenancescolors=maintenancescolors, maintenanceslinestyles=maintenanceslinestyles, plot_cens=plot_cens, plot_pb_nb_maint=plot_pb_nb_maint)
        minmaxy = min_max(minmaxy, bounds[2])
        maxx = max(maxx, (bounds[1])[2])
    end
    if ismissing(xlims)
        xlims = [0, maxx]
    end
    if ismissing(ylims)
        ylims = minmaxy
    end
    for i in 1:length(datas)
        Makie.xlims!(axs[i], xlims[1], xlims[2]*1.05)
        Makie.ylims!(axs[i], min(0,ylims[1]),ylims[2]*1.05)
    end
    
    return (xlims, ylims)
end

