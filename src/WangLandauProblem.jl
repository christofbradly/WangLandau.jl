# todo: Need possible different type for argument and output of initialise_state
# CommonSolve.jl interface
"""
    WangLandauProblem(state::StateType, moveset::MoveType)

Required user definitions:
- `initialise_state(state)::StateType`
- `histogram_size(state)::NTuple{D,Int}`
- `measure(state)::CartesianIndex`
- `random_move(state)::MoveType`
- `test_move(state, move)::CartesianIndex`
- `commit!(state, move; kwargs...)::StateType`
"""
struct WangLandauProblem{S,M}
    state::S
    moveset::M
end
# function WangLandauProblem(state, moveset)
#     return WangLandauProblem{typeof(state),typeof(moveset)}(state, moveset)
# end

"""
    CommonSolve.solve(::WangLandauProblem)::WangLandauSimulation
"""
CommonSolve.solve

"""
    initialise_state(prob::WangLandauProblem)

Initialise a new state for use in a [`WangLandauSimulation`](@ref). By
default this returns the state used to construct the
[`WangLandauProblem`](@ref), but is provided as part of the public API
with the intention of defining a lightweight struct in
`WangLandauProblem` and delaying any computationally intensive
initialisation until `WangLandauSimulation`.
"""
initialise_state(prob::WangLandauProblem) = prob.state

function histogram_size end

function measure end # -> CartesianIndex, or at least needs to index logdos

function random_move end

function test_move end

function commit! end
