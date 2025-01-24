# todo: Need possible different type for argument and output of initialise_state
# CommonSolve.jl interface
"""
    WangLandauProblem(statedef::D)

`WangLandau.jl` works by requiring user definitions of the following
generic functions:
- [`histogram_size`](@ref): the dimensions of the histogram.
- [`system_size`](@ref): the size of the system.
- [`random_trial!`](@ref): Obtain a trial move and histogram indices for
  old and new state.
- [`commit_trial!`](@ref): Upon acceptance of the trial move, update the
  state, optionally according to the indices.

Optional user definitions:
- [`initialise_state`](@ref): Optional initialisation step.
"""
struct WangLandauProblem{D}
    statedef::D
end

"""
    CommonSolve.solve(::WangLandauProblem)::WangLandauSimulation
"""
CommonSolve.solve

"""
    histogram_size(statedef::D)

This should return an integer `n`, the number of possible energy
levels accessible by `statedef`.
"""
function histogram_size end

"""
    system_size(statedef::D)

Get the canonical size of `statedef`, e.g. the number of lattice sites.
Used to determine the size of a Monte Carlo sweep.
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
initialise_state(statedef) = copy(statedef)
