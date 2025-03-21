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
- [`revert_trial!`](@ref): Upon rejection of the trial move, update the
  state, optionally according to the indices. By default returns `state`
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

Calculate a random `trial` move for `state`, as well as the `old_index`
and `new_index`.

See also [`commit_trial!`](@ref), [`revert_trial!`](@ref).
"""
function random_trial! end

"""
    commit_trial!(state::S, trial::T, old_index::I, new_index::I)

Update `state` by applying the `trial` move, using information
from `old_index` and `new_index`, if necessary. 

See also [`random_trial!`](@ref), [`revert_trial!`](@ref).
"""
function commit_trial! end

"""
    revert_trial!(state::S, trial::T, old_index::I, new_index::I)

Returns `state` to before `trial` move was performed, using information
from `old_index` and `new_index`, if necessary. 
    
By default this returns `state` unaltered. Alternatively, if this method
is defined, then it should be sufficient to define amethod for
`commit_trial!` that returns `state` unaltered.

See also [`random_trial!`](@ref), [`commit_trial!`](@ref).
"""
revert_trial!(state, trial, old_index, new_index) = state

"""
    initialise_state(statedef::DefType)::StateType

Initialise a new `state::StateType` for use in a
[`WangLandauSimulation`](@ref) based on the definition
`statedef::DefType` provided to [`WangLandauProblem`](@ref). `DefType`
and `StateType` are defined by the user and can be the same, but can be
different if initialisation is  computationally intensive. If the types
are the same then this function can simply `copy` the input, but ideally
it should reseed a new configuration. Called from
[`CommonSolve.solve!`](@ref) to seed multiple threads.
"""
function initialise_state end
