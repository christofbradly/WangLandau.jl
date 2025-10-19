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
    catchup_value(strat)
"""
catchup_value

"""
    update!(strat, sim)
"""
update!(strat::CatchupStrategy, _) = strat

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
mutable struct FixedFractionalCatchup{B} <: CatchupStrategy{B}
    const fraction::Float64
    minval::Float64
end
function FixedFractionalCatchup(f)
    return FixedFractionalCatchup{true}(f, 0.0)
end

function catchup_value(strat::FixedFractionalCatchup)
    return strat.minval * strat.fraction
end

function update!(strat::FixedFractionalCatchup, sim)
    strat.minval = minimum(sim.logdos[sim.logdos .> 0])
    return strat
end

# todo: How to implement this?
"""
    DynamicFractionalCatchup() <: CatchupStrategy

When a new state is visited for the first time, the density of states is
set to a fraction of the smallest current non-zero value. The fraction
is determined from the current value of `\\log f`.
"""
mutable struct DynamicFractionalCatchup{B} <: CatchupStrategy{B}
    value::Float64
end
function DynamicFractionalCatchup()
    return DynamicFractionalCatchup{true}(0.0)
end

function catchup_value(strat::DynamicFractionalCatchup)
    return strat.value
end

function update!(strat::DynamicFractionalCatchup; sim)
    minval = minimum(sim.logdos[sim.logdos .> 0])
    logf = current_value(sim.logf_strategy)
    strat.value = minval * exp(logf)
    return strat
end
