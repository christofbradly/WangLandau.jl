using WangLandau
using Test
using Random
using CommonSolve

include(joinpath(@__DIR__, "..", "examples", "ising.jl"))

@testset "WangLandauSimulation" begin
    Random.seed!(12345)

    L = 5
    periodic = false
    statedefn = Ising2D(L; periodic)

    sim = WangLandau.WangLandauSimulation(statedefn)
    @test isa(sim, WangLandau.WangLandauSimulation)
    @test isa(sim.logf_strategy, WangLandau.ReduceByFactor)
    @test isa(sim.flat_strategy, WangLandau.FractionOfMean)           
    @test isa(sim.catchup_strategy, WangLandau.NoCatchup)    

    sim_steps = WangLandau.WangLandauSimulation(statedefn; max_total_steps = 1)
    @test sim_steps.max_total_steps > 1
       
    sim_logf = WangLandau.WangLandauSimulation(statedefn; final_logf = 1e-3)
    @test WangLandau.final_value(sim_logf.logf_strategy) == 1e-3

    fs = WangLandau.FractionOfMean(0.75)
    sim_flat = WangLandau.WangLandauSimulation(statedefn; flat_strategy = fs)
    @test sim_flat.flat_strategy === fs

    @test_throws ArgumentError WangLandau.WangLandauSimulation(statedefn; tasks_per_thread = -1)

    io = IOBuffer()
    show(io, sim)
    out = String(take!(io))
    @test occursin("WangLandauSimulation", out)
    @test occursin("log(f)", out)
    @test occursin("iterations", out)
    @test occursin("total steps", out)
    @test occursin("elapsed time", out)
end    

@testset "CommonSolve" begin
    using Random
    Random.seed!(1234)

    L = 5
    statedefn = Ising2D(L; periodic = false)
    prob = WangLandauProblem(statedefn)

    sim = CommonSolve.init(prob; check_sweeps = 10, final_logf = 1e-3)
    @test isa(sim, WangLandau.WangLandauSimulation)
    @test sim.check_steps == 10 * WangLandau.system_size(statedefn)

    state, old_index = initialise_state(statedefn)
    logdos = zeros(Float64, WangLandau.histogram_size(statedefn))
    histogram = zeros(Int, WangLandau.histogram_size(statedefn))
    logf = 0.1
    catchup = WangLandau.NoCatchup()

    new_index = WangLandau.wl_trial!(state, old_index, statedefn, logdos, histogram, logf, catchup)
    @test 1 ≤ new_index ≤ length(histogram)
    @test any(histogram .> 0)  
    @test any(logdos .> 0)     

    hist = zeros(Int, size(sim.samples))
    task_samples = zeros(Int, max(1, sim.tasks_per_thread * Threads.nthreads()))
    chunk_size = max(1, sim.check_steps ÷ length(task_samples))

    flag = CommonSolve.step!(sim, hist, task_samples, chunk_size)
    @test isa(flag, Bool)
    @test sim.flat_checks ≥ 1
    @test sim.total_steps ≥ sim.check_steps

    sim_short = CommonSolve.init(prob; check_sweeps = 5, final_logf = 1e-2, max_total_steps = 200)
    sim_done = CommonSolve.solve!(sim_short)
    @test isa(sim_done, WangLandau.WangLandau.WangLandauSimulation)
    @test sim_done.total_steps > 0
    @test sim_done.elapsed_time ≥ 0.0
end

@testset "CatchupStrategy" begin
    using Random
    Random.seed!(1234)

    L = 4
    statedefn = Ising2D(L; periodic=false)

    # FixedFractionalCatchup 
    ffc = WangLandau.FixedFractionalCatchup(0.5)
    @test WangLandau.catchup_enabled(ffc)
    sim = WangLandau.WangLandauSimulation(statedefn)

    sim.logdos .= [0.0, 1.0, 2.0]
    WangLandau.update!(ffc, sim)
    @test ffc.minval == 1.0
    @test WangLandau.catchup_value(ffc) == 0.5
    @test ffc isa WangLandau.FixedFractionalCatchup{true}

    # DynamicFractionalCatchup
    dfc = WangLandau.DynamicFractionalCatchup()
    @test WangLandau.catchup_enabled(dfc)

    sim2 = WangLandau.WangLandauSimulation(statedefn)
    sim2.logdos .= [0.0, 2.0, 3.0]
    WangLandau.update!(dfc, sim2)
    @test dfc.value ≥ 0.0
    @test WangLandau.catchup_value(dfc) == dfc.value

    # FFC trial
    state, idx = initialise_state(statedefn)
    logdos = zeros(Float64, WangLandau.histogram_size(statedefn))
    hist = zeros(Int, WangLandau.histogram_size(statedefn))
    logf = 0.1

    new_idx = WangLandau.wl_trial!(state, idx, statedefn, logdos, hist, logf, ffc)
    @test 1 ≤ new_idx ≤ length(hist)
    @test any(logdos .> 0)
    @test any(hist .> 0)

    # DFC trial
    state, idx = initialise_state(statedefn)
    logdos = zeros(Float64, WangLandau.histogram_size(statedefn))
    hist2 = zeros(Int, WangLandau.histogram_size(statedefn))
    logf = 0.1

    new_idx = WangLandau.wl_trial!(state, idx, statedefn, logdos, hist2, logf, dfc)
    @test 1 ≤ new_idx ≤ length(hist)
    @test any(logdos .> 0)
    @test any(hist .> 0)

    sim_fixed = WangLandau.WangLandauSimulation(statedefn; catchup_strategy = ffc)
    @test sim_fixed.catchup_strategy === ffc

    sim_dynamic = WangLandau.WangLandauSimulation(statedefn; catchup_strategy = dfc)
    @test sim_dynamic.catchup_strategy === dfc

    hist3 = zeros(Int, size(sim_fixed.samples))
    task_samples1 = zeros(Int, max(1, sim_fixed.tasks_per_thread * Threads.nthreads()))
    chunk_size1 = max(1, sim_fixed.check_steps ÷ length(task_samples1))
    CommonSolve.step!(sim_fixed, hist3, task_samples1, chunk_size1)
    @test sim_fixed.flat_checks ≥ 1

    hist4 = zeros(Int, size(sim_dynamic.samples))
    task_samples2 = zeros(Int, max(1, sim_dynamic.tasks_per_thread * Threads.nthreads()))
    chunk_size2 = max(1, sim_dynamic.check_steps ÷ length(task_samples2))
    CommonSolve.step!(sim_dynamic, hist4, task_samples2, chunk_size2)
    @test sim_dynamic.flat_checks ≥ 1
end