# todo: Need possible different type for argument and output of initialise_state
# CommonSolve.jl interface
"""
    WangLandauProblem(state::S)

`WangLandau.jl` works by requiring user definitions of the following
generic functions:
- [`histogram_size`](@ref): get the dimensions of the histogram.
- [`random_trial!`](@ref): Obtain a trial move and histogram indices for
  old and new state.
- [`commit_trial!`](@ref): Upon acceptance of the trial move, update the
  state, optionally according to the indices.

Optional user definitions:
- [`initialise_state`](@ref): Optional initialisation step.
"""
struct WangLandauProblem{S}
    state::S
end

"""
    CommonSolve.solve(::WangLandauProblem)::WangLandauSimulation
"""
CommonSolve.solve

"""
    histogram_size(state::S)
"""
function histogram_size end

"""
    random_trial!(state::S) -> trial::T, old_index::I, new_index::I

See also [`commit_trial!`](@ref).
"""
function random_trial! end

"""
    commit_trial!(state::S, trial::T, old_index::I, new_index::I)

See also [`random_trial!`](@ref).
"""
function commit_trial! end

"""
    initialise_state(state::S)

Initialise a new `state` for use in a [`WangLandauSimulation`](@ref). By
default this returns the state used to construct the
[`WangLandauProblem`](@ref), but is provided as part of the public API
with the intention of defining a lightweight struct in
`WangLandauProblem` and delaying any computationally intensive
initialisation until `WangLandauSimulation`. Called from
[`CommonSolve.init`](@ref).
"""
initialise_state(state) = state
