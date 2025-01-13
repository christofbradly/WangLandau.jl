# CommonSolve.jl interface
"""
    WangLandauProblem()
"""
struct WangLandauProblem{}

end
function WangLandauProblem(n; kwargs...)

    return WangLandauProblem()
end

"""
    WangLandauSimulation()
"""
struct WangLandauSimulation{S}
    state::S
end
function WangLandauSimulation(state; kwargs...)
    
    return WangLandauSimulation{typeof(state)}(state)
end

"""
    CommonSolve.init(problem::WangLandauProblem) -> WangLandauSimulation

Initialise a [`WangLandauSimulation`](@ref) based on `problem`.
"""
function CommonSolve.init(prob::WangLandauProblem)
    state = initialise(prob.state)
    return WangLandauSimulation{S}(state)
end

"""
    CommonSolve.solve!(sim::WangLandauSimulation; kwargs...)

Run WangLandau algorithm on `sim`.
"""
function CommonSolve.solve!(sim::WangLandauSimulation)

    return sim
end

"""
    initialise(state)

Initialise a `state` for use in a [`WangLandauSimulation`](@ref). By
default this is the identity function, but is provided as part of the
public API with the intention of defining a lightweight struct in
[`WangLandauProblem`](@ref) and delaying any computationally intensive
initialisation until `WangLandauSimulation`.
"""
initialise(state) = state