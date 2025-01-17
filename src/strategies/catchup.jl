"""
    CatchupStrategy{B}

Strategy for determining how to update the density of states when a new
state is visited for the first time.

This occurs *before* the density of states is updated with parameter `f`.

API: 
- [`catchup_enabled(strat::CatchupStrategy)`](@ref)
- [`catchup_value(strat::CatchupStrategy, sim)`](@ref): only called if
`B == true`
"""
abstract type CatchupStrategy{B} end

"""
    catchup_enabled(strat)
"""
catchup_enabled(::CatchupStrategy{B}) where {B} = B

"""
    catchup_value(strat, sim)
"""
catchup_value

"""
    NoCatchup() <: CatchupStrategy{false}

The default strategy, does nothing.
"""
struct NoCatchup{B} <: CatchupStrategy{B}
end
NoCatchup() = NoCatchup{false}()
function catchup_value(::NoCatchup)
    return 0.0
end

"""
    FixedFractionalCatchup(f) <: CatchupStrategy

When a new state is visited for the first time, the density of states is
set to a fixed fraction of the smallest current non-zero value.
"""
struct FixedFractionalCatchup{B} <: CatchupStrategy{B}
    fraction::Float64
end
function FixedFractionalCatchup(f)
    return FixedFractionalCatchup{true}(f)
end

function catchup_value(strat::FixedFractionalCatchup; sim)
    dos = sim.logdos
    return minimum(dos[dos .> 0]) * strat.fraction
end

# todo: How to implement this?
"""
    DynamicFractionalCatchup() <: CatchupStrategy

When a new state is visited for the first time, the density of states is
set to a fraction of the smallest current non-zero value. The fraction
is determined from the current value of `\\log f`.
"""
struct DynamicFractionalCatchup{B} <: CatchupStrategy{B}
    DynamicFractionalCatchup() = new{true}()
end

function catchup_value(::DynamicFractionalCatchup; sim)
    dos = sim.logdos
    logf = current_value(sim.logf_strategy)
    return minimum(dos[dos .> 0]) * (1 - logf)
end
