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
    update!(strat, sim)

Optionally update `strat` with information from the simulation `sim`.
Returns `nothing`.
"""
update!(::FlatHistogramStrategy, _) = nothing

"""
    FractionOfMean(tol)

Define the histogram flatness criterion to be when the smallest non-zero
value is greater than `tol` times the mean (of non-zero values).
"""
struct FractionOfMean <: FlatHistogramStrategy
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

"""
    StableNumVisits(min, checksteps)

Define the histogram flatness criterion to be when the number of visited
states remains constant for a number of steps `checksteps`. A state is
considered visited if it has been sampled at least `minvisits` times.

This strategy may be more effective for two-dimensional state spaces
(i.e. energy and order parameter), see e.g. [Tsai, Wang & Landau, *Braz.
J. Phys.* **38** 2008
](https://doi.org/10.1590/S0103-97332008000100003).

Suggested parameters are `min = 2000` and `check = N * 10^6` for
system size `N`.

"""
mutable struct StableNumVisits <: FlatHistogramStrategy
    const min::Int
    const checksteps::Int
    numvisits::Int
    lastcheck::Int
    stablesteps::Int
    iter::Int
    function StableNumVisits(min, checksteps)
        0 < tol < 1 || throw(ArgumentError("Tolerance must be less than 1."))
        return new(min, checksteps, 0, 0, 0, 0)
    end    
end

function isflat(strat::StableNumVisits, hist)
    new_numvisits = length(hist[hist .> strat.min])
    if new_numvisits > strat.numvisits
        strat.stablesteps = 0
        strat.numvisits = new_numvisits
        return false
    end
    # numvisits is stable, now check for duration
    return strat.stablesteps > checksteps
end

function update!(strat::StableNumVisits, sim)
    if sim.flat_iterations > strat.iter
        strat.iter = sim.flat_iterations
        strat.stablesteps = 0
        strat.numvisits = 0
        return strat
    end
    strat.stablesteps = sim.total_steps - strat.lastcheck
    return strat
end
