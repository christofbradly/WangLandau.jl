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

    sim = WangLandauProblem(statedefn)
    @test sim isa WangLandau.WangLandauSimulation
    @test sim.logf_strategy isa WangLandau.ReduceByFactor
    @test sim.flat_strategy isa WangLandau.FractionOfMean
    @test sim.catchup_strategy isa WangLandau.NoCatchup

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
    @test sim isa WangLandau.WangLandauSimulation

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