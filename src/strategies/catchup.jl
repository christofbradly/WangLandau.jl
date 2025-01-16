"""
    CatchupStrategy{B}

Strategy for determining how to catchup the estimate of the density of
states when a new state is visited for the first time.

API:
`catchup(strat::CatchupStrategy; kwargs...)`: only called if `B == true`
"""
abstract type CatchupStrategy{B} end

"""
    NoCatchup() <: CatchupStrategy{false}

The default strategy, does nothing.
"""
struct NoCatchup <: CatchupStrategy{false}
end
function catchup(::NoCatchup)
    return zero(T)
end

"""
    FixedFractionalCatchup(f) <: CatchupStrategy

When a new state is visited for the first time, the density of states is
set to a fixed fraction of the smallest current non-zero value.
"""
struct FixedFractionalCatchup{B} <: CatchupStrategy
    fraction::Float64
end
function FixedFractionalCatchup(f)
    return FixedFractionalCatchup{true}(f)
end

function catchup(strat::FixedFractionalCatchup; dos_data)
    return minimum(dos_data[dos_data .> 0]) * strat.fraction
end

# todo: How to implement this?
"""
    DynamicFractionalCatchup() <: CatchupStrategy

When a new state is visited for the first time, the density of states is
set to a fraction of the smallest current non-zero value. The fraction
is determined from the current value of `\\log f`.
"""
struct DynamicFractionalCatchup{B} <: CatchupStrategy
    fraction::Float64
end
function DynamicFractionalCatchup()
    return DynamicFractionalCatchup{true}()
end

function catchup(strat::DynamicFractionalCatchup; dos_data, state)
    return minimum(dos_data[dos_data .> 0]) * (1 - logf)
end
