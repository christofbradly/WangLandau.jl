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
value is greater than `tol` times the mean, which is calculated from
states that have at least `min` visits.
"""
@kwdef struct FractionOfMean <: FlatHistogramStrategy
    tol::Float64
    min::Int = 0
    function FractionOfMean(tol, min)
        0 < tol < 1 || throw(ArgumentError("Tolerance must be less than 1."))
        min ≥ 0 || throw(ArgumentError("Min must be non-negative."))
        return new(tol, min)
    end
end
FractionOfMean(tol; kwargs...) = FractionOfMean(; tol, kwargs...)

function isflat(strat::FractionOfMean, hist)
    min = strat.min
    numnonzeros = length(hist[hist .> min])
    # avoid false positive during initial warmup
    numnonzeros == 0 && throw(ErrorException("No samples in local histogram."))
    numnonzeros == 1 && return false
    mean = sum(hist[hist .> min]) / numnonzeros
    gap = minimum(hist[hist .> min]) / mean
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
