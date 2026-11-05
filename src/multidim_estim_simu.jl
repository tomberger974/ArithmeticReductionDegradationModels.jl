function tests(true_mvw::MvWienerAR, degradationdata::DegradationData, N::Int64)
    estimations = Vector{MvWienerAR}(undef, N)

    t1 = time()

    Threads.@threads for i in 1:N
        degradationdata_i = rand(true_mvw, degradationdata)

        estimations[i] = fit_mle(true_mvw, degradationdata_i)
    end

    println("Elapsed time: ", time() - t1)

    return estimations
end

function param(estimation_set::ArithmeticReductionDegradationModels.MvWienerAR,prop::Symbol;sel::Union{Symbol, Tuple{Union{String, Symbol}, Union{String, Symbol}}})
    return get(getproperty(estimation_set,prop),sel,nothing)
end