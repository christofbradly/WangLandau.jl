"""
    WangLandauSimulation()

Keyword arguments:
- `check_steps = 1000`:
- `max_total_steps`:
- `final_logf = 1e-6`: when ``\\log f`` reaches this value the
  simulation ends.
- `logf_strategy = LogReduceByFactor(; final = final_logf)`: Controls
  how ``f`` is updated. Overrides `final_logf`.
- `tol = 0.8`: Set tolerance for the flatness of the histogram.
- `flat_strategy = FractionOfMean(tol)`: Define the flatness criterion
  for the histogram. Overrides `tol`.
"""
mutable struct WangLandauSimulation{S,M,D,C}
    state::S
    moveset::M
    samples::Array{Int,D}
    logdos::Array{Float64,D}
    logf_strategy::DosIncrementStrategy
    flat_strategy::FlatHistogramStrategy
    check_steps::Int
    flat_checks::Int
    flat_iterations::Int
    total_steps::Int
    max_total_steps::Int
    catchup::CatchupStrategy{C}
    elapsed_time::Float64
end
function WangLandauSimulation(state::S,moveset::M;
    check_steps = 1000,
    max_total_steps = 1e6,
    final_logf = 1e-6,
    logf_strategy = nothing,
    flat_tolerance = 0.9,
    flat_strategy = nothing,
    catchup = NoCatchup()
    ) where {S,M}

    dims = histogram_size(prob.state)
    logdos = zeros(dims)
    samples = zeros(Int, dims)
    D = length(dims)
    C = catchup_enabled(catchup)

    if isnothing(logf_strategy)
        logf_strategy = LogReduceByFactor(final_logf)
    end
    if isnothing(flat_strategy)
        flat_strategy = FractionOfMean(tol)
    end
    if max_total_steps < expected_iterations(flat_strategy) * check_steps
        @warn "max_total_steps is less than expected number; overwriting to expected number of steps."
        max_total_steps = expected_iterations(flat_strategy) * check_steps
    end

    return WangLandauSimulation{S,M,D,C}(
        state,
        moveset,
        samples,
        logdos,
        logf_strategy,
        flat_strategy,
        check_steps,
        0,
        0,
        0,
        max_total_steps,
        catchup,
        0.0
        )
end

function Base.show(io::IO, sim::WangLandauSimulation{S}) where {S}
    print(io, "WangLandauSimulation{", S, "}")
    print(io, "log(f) = ", sim.logf, " (final: ", sim.final_logf, ")")
    print(io, "iterations: ", sim.flat_iterations, "(checks: ", sim.flat_checks,")")
    print(io, "total_steps: ", sim.total_steps)
    return nothing
end

"""
    CommonSolve.init(problem::WangLandauProblem) -> WangLandauSimulation

Initialise a [`WangLandauSimulation`](@ref) based on `problem`.
"""
function CommonSolve.init(prob::WangLandauProblem)
    return WangLandauSimulation(prob.state, prob.moveset; kwargs...)
end

"""
    CommonSolve.step!(sim::WangLandauSimulation, temp_hist)

Run `sim` for a single step, using `temp_hist` as the local histogram.
"""
function CommonSolve.step!(
    sim::WangLandauSimulation{<:Any,<:Any,<:Any,C},
    temp_hist,
    ) where {C}

    (; state, flat_iterations, logdos, logf_strategy) = sim

    old_index = measure(state)
    move = random_move(state)
    new_index = test_move(state, move)
    
    old_dos = logdos[old_index]
    new_dos = logdos[new_index]

    if log(rand()) < old_dos - new_dos
        commit!(state, move, new_index)
    else
        new_index = old_index
    end
    temp_hist[new_index] += 1

    if C && iszero(logdos[new_index]) && flat_iterations > 0
        logdos[new_index] = catchup_value(sim.catchup; logdos)
    else
        logdos[new_index] += current(logf_strategy)
    end

    if iszero(logdos[new_index])
        if C && flat_iters > 0
            logdos[new_index] = catchup_value(catchup; logdos)
        else
            logdos[new_index] += current(logf_strategy)
        end
    else
        logdos[new_index] += current(logf_strategy)
    end
    return nothing
end

"""
    CommonSolve.solve!(sim::WangLandauSimulation; kwargs...)

Run WangLandau algorithm on `sim`.
"""
function CommonSolve.solve!(sim::WangLandauSimulation)
    state = initialise_state(sim.state)
    
    temp_hist = zeros(Int, size(sim.samples))

    num_iterations = expected_iterations(sim.logf_strategy)

    starting_time = time() + sim.elapsed_time
    @withprogress name = "WangLandau" begin
        while !isconverged(sim.logf_strategy)
            for _ in 1:sim.check_steps
                step!(sim, temp_hist)
            end

            if isflat(sim.flat_strategy, temp_hist)
                update!(sim.logf_strategy)
                sim.flat_iterations += 1
                samples .+= temp_hist
                temp_hist = zeros(Int, dims)
            end
            update!(sim.flat_strategy, sim)

            sim.flat_checks += 1
            sim.total_steps += check_steps

            @logprogress flat_iterations / num_iterations
        end
    end

    sim.elapsed_time = time() - starting_time
    return sim
end
