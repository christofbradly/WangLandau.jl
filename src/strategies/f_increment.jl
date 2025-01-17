# to-do: the current value of logf should be local to solve!, not held
# in these strategies. Then 
"""
    DosIncrementStrategy
"""
abstract type DosIncrementStrategy end

# API:
"""
    current_value(strat)

Return the current value of ``\\log f`` from the `strat`.
"""
current_value

"""
    update!(strat)    
"""
update!

"""
    isconverged(strat)
"""
isconverged

"""
    expected_iterations(strat)
"""
expected_iterations

"""
    ReduceByFactor(; kwargs...)

A [`DosIncrementStrategy`](@ref) to control parameter ``\\log f`` by
multiplying by a constant factor every time a flat histogram occurs,
until a desired minimum value is achieved.

Keyword arguments and their default values are:
- `initial = 1.0`
- `factor = 0.5`: must be less than 1
- `final = 1e-6`
"""
mutable struct ReduceByFactor <: DosIncrementStrategy
    const initial::Float64
    current::Float64
    const factor::Float64
    const final::Float64
end
function ReduceByFactor(; initial = 1.0, factor = 0.5, final = 1e-6)
    0 < factor < 1 || throw(ArgumentError("Reduction factor must be less than 1."))
    initial > final || throw(ArgumentError("Initial value must be greater than final value."))
    return ReduceByFactor(initial, initial, factor, final)    
end

current_value(strat::ReduceByFactor) = strat.current

function update!(strat::ReduceByFactor)
    strat.current *= strat.factor
    return strat
end

isconverged(strat::ReduceByFactor) = strat.current < strat.final

function expected_iterations(strat::ReduceByFactor)
    (; initial, factor, final) = strat
    n = ceil(log(final / initial) / log(factor))
    return Int(n)
end
