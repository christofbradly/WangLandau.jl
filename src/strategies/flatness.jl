"""
    FlatHistogramStrategy
"""
abstract type FlatHistogramStrategy end

"""
    isflat(strat, hist)::Bool

Checks if histogram `hist` is flat according to the
[`FlatHistogramStrategy`](@ref).
"""
isflat

"""
    FractionOfMean(tol)

Define the histogram flatness criterion to be when the smallest non-zero
value is greater than `tol` times the mean (of non-zero values).
"""
struct FractionOfMean
    tol::Float64
    function FractionOfMean(tol)
        0 < tol < 1 || throw(ArgumentError("Tolerance must be less than 1."))
        return new(tol)
    end    
end

function isflat(strat::FractionOfMean, hist)
    numnonzeros = length(hist[hist .> 0])
    # avoid false positive during initial warmup
    numnonzeros == 0 && throw(ErrorException("No samples in local histogram."))
    numnonzeros == 1 && return false
    mean = sum(hist[hist .> 0]) / numnonzeros
    gap = minimum(hist[hist .> 0]) / mean
    return gap > strat.tol
end
