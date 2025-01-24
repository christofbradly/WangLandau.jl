"""
    WangLandauSimulation()

Keyword arguments:
- `check_sweeps = 100`: The number of sweeps to perform before checking
  for flatness. A sweep is `N` steps where `N` is the size of the
  system. See [`system_size`](@ref).
- `max_total_steps = Inf`:
- `final_logf = 1e-6`: when ``\\log f`` reaches this value the
  simulation ends.
- `logf_strategy = ReduceByFactor(; final = final_logf)`: Controls how
  ``f`` is updated. Overrides `final_logf`.
- `tol = 0.8`: Set tolerance for the flatness of the histogram.
- `flat_strategy = FractionOfMean(tol)`: Define the flatness criterion
  for the histogram. Overrides `tol`.
- `tasks_per_thread = 4`: a multiplier for determining the number of
  tasks to `@spawn`. Set to `0` to disable concurrent threads even if more
  than one thread is running.
"""
mutable struct WangLandauSimulation{S,D,I,F,C}
    statedef::S
    logf_strategy::I
    flat_strategy::F
    catchup_strategy::C
    const check_steps::Int
    flat_checks::Int
    flat_iterations::Int
    total_steps::Int
    const max_total_steps::Float64
    const tasks_per_thread::Int
    elapsed_time::Float64
    samples::Array{Int,D}
    logdos::Array{Float64,D}
end
function WangLandauSimulation(statedef::S;
    check_sweeps = 100,
    max_total_steps = Inf,
    tasks_per_thread = 4,
    final_logf = 1e-6,
    logf_strategy = nothing,
    flat_tolerance = 0.9,
    flat_strategy = nothing,
    catchup_strategy = NoCatchup()
    ) where {S}

    dims = histogram_size(statedef)
    logdos = zeros(dims)
    samples = zeros(Int, dims)
    D = length(dims)
    # C = catchup_enabled(catchup_strategy)

    check_steps = check_sweeps * system_size(statedef)

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
    tasks_per_thread ≥ 0 || throw(ArgumentError("Task multiplier must be non-negative."))
    I, F, C = typeof.((logf_strategy, flat_strategy, catchup_strategy))

    return WangLandauSimulation{S,D,I,F,C}(
        statedef,
        logf_strategy,
        flat_strategy,
        catchup_strategy,
        check_steps,
        0,
        0,
        0,
        max_total_steps,
        tasks_per_thread,
        0.0,
        samples,
        logdos,
        )
end

function Base.show(io::IO, sim::WangLandauSimulation)
    logf = current_value(sim.logf_strategy)
    final_logf = final_value(sim.logf_strategy)
    println(io, "WangLandauSimulation(", sim.statedef, ")")
    println(io, "  log(f) = ", logf, " (final: ", final_logf, ")")
    println(io, "  iterations: ", sim.flat_iterations, " (checks: ", sim.flat_checks,")")
    println(io, "  total steps: ", sim.total_steps)
    @printf io "  elapsed time: %.2e s" sim.elapsed_time
    # print(io, "  elapsed time: ", sim.elapsed_time, " s")
    return nothing
end

"""
    CommonSolve.init(problem::WangLandauProblem; kwargs...) -> WangLandauSimulation

Initialise a [`WangLandauSimulation`](@ref) based on `problem`.
"""
function CommonSolve.init(prob::WangLandauProblem; kwargs...)
    return WangLandauSimulation(prob.statedef; kwargs...)
end

"""
    wl_trial!(state, logdos, temp_hist, logf, catchup)

Obtain a single trial move, compare to current `state` and commit or
reject Then increment the density of states `logdos` and histogram
`temp_hist`, with `logf` and `1`, respectively.
"""
function wl_trial!(state, logdos, histogram, logf, catchup_strategy::CatchupStrategy{C}) where {C}

    trial, old_index, new_index = random_trial!(state)
    
    old_dos = Atomix.@atomic logdos[old_index]
    new_dos = Atomix.@atomic logdos[new_index]

    # faster than log(rand())
    if rand() < exp(old_dos - new_dos)
        commit_trial!(state, trial, old_index, new_index)
    else
        new_index = old_index
    end
    
    if C && iszero(new_dos)
        logdos_incr = catchup_value(catchup_strategy)   # not atomic
    else
        logdos_incr = logf
    end
    Atomix.@atomic logdos[new_index] += logdos_incr
    Atomix.@atomic histogram[new_index] += 1

    return nothing
end

"""
    CommonSolve.step!(sim::WangLandauSimulation, histogram)

Run `sim` for a single iteration until the `histogram` is flat.
"""
function CommonSolve.step!(sim::WangLandauSimulation, histogram, states, task_samples, chunk_size)
    (; logdos, logf_strategy, catchup_strategy) = sim
    logf = current_value(logf_strategy)

    chunks = Base.Iterators.partition(1:sim.check_steps, chunk_size)
    tasks = map(enumerate(chunks)) do (i, chunk)
        Threads.@spawn begin
            state = states[$i]
            for _ in chunk
                wl_trial!(state, logdos, histogram, logf, catchup_strategy)
            end
            return length(chunk)
        end
    end
    task_samples .+= fetch.(tasks)

    flat = isflat(sim.flat_strategy, histogram)
    if flat
        update!(sim.logf_strategy)
        sim.flat_iterations += 1
        sim.samples .+= histogram
        histogram .= 0
        @info "flat!", sim.flat_iterations, logf, sim.flat_checks, task_samples
    end
    update!(sim.flat_strategy, sim)
    update!(sim.catchup_strategy, sim)
    sim.flat_checks += 1
    sim.total_steps += sim.check_steps

    return flat
end

"""
    CommonSolve.solve!(sim::WangLandauSimulation; kwargs...)

Run WangLandau algorithm on `sim`.
"""
function CommonSolve.solve!(sim::WangLandauSimulation)
    # multi-threading setup
    task_sets = max(1, sim.tasks_per_thread * Threads.nthreads())
    chunk_size = sim.check_steps ÷ task_sets + 1
    states = [initialise_state(sim.statedef) for _ in 1:task_sets]
    task_samples = zeros(Int, task_sets)

    # local diagnostics
    temp_hist = zeros(Int, size(sim.samples))
    total_iterations = expected_iterations(sim.logf_strategy)

    starting_time = time() + sim.elapsed_time
    @info "Starting simulation..."
    @withprogress name = "WangLandau" begin
        while !isconverged(sim.logf_strategy)
            flag = CommonSolve.step!(sim, temp_hist, states, task_samples, chunk_size)
            if flag
                @debug "Flat! after", sim.flat_iterations, " iterations."
            end

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
