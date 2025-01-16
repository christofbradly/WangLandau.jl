abstract type DosIncrementStrategy end

# API:
# f_value
# update!
# isconverged
# expected

"""
    LogReduceByFactor(; kwargs...)

A [`DosIncrementStrategy`](@ref) to control parameter ``\\log f`` by
multiplying by a constant factor every time a flat histogram occurs,
until a desired minimum value is achieved.

Keyword arguments and their default values are:
- `initial = 1.0`
- `factor = 0.5`: must be less than 1
- `final = 1e-6`
"""
mutable struct LogReduceByFactor <: DosIncrementStrategy
    const initial::Float64
    current::Float64
    const factor::Float64
    const final::Float64
end
function LogReduceByFactor(; initial = 1.0, factor = 0.5, final = 1e-6)
    0 < factor < 1 || throw(ArgumentError("Reduction factor must be less than 1."))
    initial > final || throw(ArgumentError("Initial value must be greater than final value."))
    return LogReduceByFactor(initial, initial, factor, final)    
end

f_value(strat::LogReduceByFactor) = st.current

function update!(strat::LogReduceByFactor)
    strat.current *= strat.factor
    return strat
end

isconverged(strat::LogReduceByFactor) = strat.current < strat.final

function expected_iterations(strat::LogReduceByFactor)
    (; initial, factor, final) = strat
    n = ceil(log(final / initial) / log(factor))
    return Int(n)
end
