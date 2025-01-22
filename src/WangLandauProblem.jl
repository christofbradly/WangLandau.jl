# todo: Need possible different type for argument and output of initialise_state
# CommonSolve.jl interface
"""
    WangLandauProblem(statedef::S)

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
    statedef::S
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
    system_size(state::S)
"""
function system_size end

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
    initialise_state(statedef::DefType)::StateType

Initialise a new `state::StateType` for use in a
[`WangLandauSimulation`](@ref) based on the defintion
`statedef::DefType` provided to [`WangLandauProblem`](@ref). `DefType`
and `StateType` are defined by the user and can be the same. By default
this function copies `statedef`, but is provided as part of the public
API with the intention of defining a lightweight struct in
`WangLandauProblem` and delaying any computationally intensive
initialisation until [`solve`](@ref) is called. It also allows for
reseeding the state for use with multiple threads. Called from
[`CommonSolve.init`](@ref).
"""
initialise_state(state) = copy(state)
