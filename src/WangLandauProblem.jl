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
struct WangLandauSimulation{}

end
function WangLandauSimulation(prob::WangLandauProblem; kwargs...)
    return WangLandauSimulation()
end

"""
    CommonSolve.init(problem::WangLandauProblem) -> WangLandauSimulation

Initialise a [`WangLandauSimulation`](@ref) based on `problem`.
"""
function CommonSolve.init(prob::WangLandauProblem)

    return WangLandauSimulation(prob)
end

"""
    CommonSolve.solve!(sim::WangLandauSimulation; kwargs...)

Run WangLandau algorithm on `sim`.
"""
function CommonSolve.solve!(sim::WangLandauSimulation)

    return sim
end
