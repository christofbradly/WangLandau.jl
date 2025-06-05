# todo: Need possible different type for argument and output of initialise_state
# CommonSolve.jl interface
"""
    WangLandauProblem(statedefn::D)

This is a simple wrapper for the user-defined `statedefn` that enables
the `CommonSolve.jl` interface. It can be overloaded so that `statedefn`
does not have to be instantiated directly first.
"""
struct WangLandauProblem{D}
    statedefn::D
end

"""
    CommonSolve.solve(::WangLandauProblem)::WangLandauSimulation
"""
CommonSolve.solve

"""
    histogram_size(statedefn::D)

Return a `Tuple` of integers that specify the size of the histograms.
The first is canonically the number of possible energy levels accessible
by `statedefn`.
"""
function histogram_size end

"""
    system_size(statedefn::D)

Get the canonical size of `statedefn`, e.g. the number of lattice sites.
Used to determine the size of a Monte Carlo sweep.
"""
function system_size end

"""
    random_trial!(state::S, statedefn::D) -> trial::T

Calculate a random `trial` move for `state` based on `statedefn`.

If a trial move cannot be found returns `nothing`.

See also [`commit_trial!`](@ref), [`revert_trial!`](@ref).
"""
function random_trial! end

"""
    hist_index(state::S, statedefn::D, trial::T, old_index::I) -> new_index::I

Calculate the `new_index` for accessing the density of states, using
`statedefn`, `trial` or `old_index`, if necessary.

See also [`random_trial`](@ref).
"""
function hist_index end

"""
    commit_trial!(state::S, statedefn::D, trial::T, old_index::I, new_index::I)

Update `state` by applying the `trial` move, using information
from `old_index` and `new_index`, if necessary. 

See also [`random_trial!`](@ref), [`revert_trial!`](@ref).
"""
function commit_trial! end

"""
    revert_trial!(state::S, statedefn::D, trial::T, old_index::I, new_index::I)

Returns `state` to before `trial` move was performed, using information
from `old_index` and `new_index`, if necessary. 
    
By default this returns `state` unaltered. Alternatively, if this method
is defined, then it should be sufficient to define amethod for
`commit_trial!` that returns `state` unaltered.

See also [`random_trial!`](@ref), [`commit_trial!`](@ref).
"""
revert_trial!(state, _, _, _, _) = state

"""
    initialise_state(statedefn::D) -> (state::S, index::I)

Initialise a new `state::S` for use in a [`WangLandauSimulation`](@ref)
based on the definition `statedefn::D` provided to
[`WangLandauProblem`](@ref). `D` and `S` are defined by the user and `D`
should be immutable. This function could copy a configuration that is
stored in `statedefn`, but ideally it should reseed a new configuration.
Called from [`CommonSolve.solve!`](@ref) to seed multiple threads.
"""
function initialise_state end
