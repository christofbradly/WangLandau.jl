"""
    WangLandauSimulation()

Keyword arguments:
- `check_steps = 1000`:
- `max_total_steps = Inf`:
- `final_logf = 1e-6`: when ``\\log f`` reaches this value the
  simulation ends.
- `logf_strategy = ReduceByFactor(; final = final_logf)`: Controls
  how ``f`` is updated. Overrides `final_logf`.
- `tol = 0.8`: Set tolerance for the flatness of the histogram.
- `flat_strategy = FractionOfMean(tol)`: Define the flatness criterion
  for the histogram. Overrides `tol`.
"""
mutable struct WangLandauSimulation{S,C,D}
    state::S
    samples::Array{Int,D}
    logdos::Array{Float64,D}
    logf_strategy::DosIncrementStrategy
    flat_strategy::FlatHistogramStrategy
    check_steps::Int
    flat_checks::Int
    flat_iterations::Int
    total_steps::Int
    max_total_steps::Float64
    catchup::CatchupStrategy{C}
    elapsed_time::Float64
end
function WangLandauSimulation(state::S;
    check_steps = 1000,
    max_total_steps = Inf,
    final_logf = 1e-6,
    logf_strategy = nothing,
    flat_tolerance = 0.9,
    flat_strategy = nothing,
    catchup = NoCatchup()
    ) where {S}

    dims = histogram_size(state)
    logdos = zeros(dims)
    samples = zeros(Int, dims)
    D = length(dims)
    C = catchup_enabled(catchup)

    if isnothing(logf_strategy)
        logf_strategy = ReduceByFactor(; final = final_logf)
    end
    if isnothing(flat_strategy)
        flat_strategy = FractionOfMean(flat_tolerance)
    end
    if max_total_steps < expected_iterations(logf_strategy) * check_steps
        @warn "`max_total_steps` is less than expected number; overwriting to expected number of steps."
        max_total_steps = expected_iterations(logf_strategy) * check_steps
    end

    return WangLandauSimulation{S,C,D}(
        state,
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
    logf = current_value(sim.logf_strategy)
    final_logf = final_value(sim.logf_strategy)
    println(io, "WangLandauSimulation{", S, "}")
    println(io, "  log(f) = ", logf, " (final: ", final_logf, ")")
    println(io, "  iterations: ", sim.flat_iterations, " (checks: ", sim.flat_checks,")")
    print(io, "  total steps: ", sim.total_steps)
    return nothing
end

"""
    CommonSolve.init(problem::WangLandauProblem; kwargs...) -> WangLandauSimulation

Initialise a [`WangLandauSimulation`](@ref) based on `problem`.
"""
function CommonSolve.init(prob::WangLandauProblem; kwargs...)
    state = initialise_state(prob.state)
    return WangLandauSimulation(state; kwargs...)
end

"""
    CommonSolve.step!(sim::WangLandauSimulation, temp_hist)

Run `sim` for a single step, using `temp_hist` as the local histogram.
"""
function CommonSolve.step!(
    sim::WangLandauSimulation{<:Any,C},
    temp_hist,
    ) where {C}

    (; state, flat_iterations, logdos, logf_strategy) = sim

    trial, old_index, new_index = random_trial!(state)
    
    old_dos = logdos[old_index]
    new_dos = logdos[new_index]

    # faster than log(rand())
    if rand() < exp(old_dos - new_dos)
        commit_trial!(state, trial, old_index, new_index)
    else
        new_index = old_index
    end
    temp_hist[new_index] += 1

    if C && iszero(logdos[new_index]) && flat_iterations > 0
        logdos[new_index] = catchup_value(sim.catchup; logdos)
    end
    logdos[new_index] += current_value(logf_strategy)

    return nothing
end

"""
    CommonSolve.solve!(sim::WangLandauSimulation; kwargs...)

Run WangLandau algorithm on `sim`.
"""
function CommonSolve.solve!(sim::WangLandauSimulation)
    dims = size(sim.samples)
    temp_hist = zeros(Int, dims)

    total_iterations = expected_iterations(sim.logf_strategy)

    starting_time = time() + sim.elapsed_time
    @info "Starting simulation..."
    @withprogress name = "WangLandau" begin
        while !isconverged(sim.logf_strategy)
            for _ in 1:sim.check_steps
                step!(sim, temp_hist)
            end

            if isflat(sim.flat_strategy, temp_hist)
                update!(sim.logf_strategy)
                sim.flat_iterations += 1
                sim.samples .+= temp_hist
                temp_hist = zeros(Int, dims)
            end
            update!(sim.flat_strategy, sim)

            sim.flat_checks += 1
            sim.total_steps += sim.check_steps

            if sim.total_steps > sim.max_total_steps 
                @warn "Global maximum number of steps reached, ending simulation."
                break
            end
            
            @logprogress sim.flat_iterations / total_iterations
        end
    end

    sim.elapsed_time = time() - starting_time
    @info "... done!"
    return sim
end
