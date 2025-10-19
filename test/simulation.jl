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

@testset "FixedFractionalCatchup" begin
    Random.seed!(12345)

    L = 4
    statedefn = Ising2D(L; periodic=false)
    prob = WangLandauProblem(statedefn)

    ffc = WangLandau.FixedFractionalCatchup(0.25)
    sim1 = CommonSolve.init(prob; catchup_strategy = ffc)
    sim1.logdos .= [1e-8, 2e-8, 3e-8, 4e-8][mod1.(1:length(sim1.logdos), 4)]

    oldmin = ffc.minval
    WangLandau.update!(ffc, sim1)
    @test ffc.minval != oldmin
    @test WangLandau.catchup_enabled(ffc) == true
    @test WangLandau.catchup_value(ffc) ≈ ffc.minval * ffc.fraction

    state, old_index = initialise_state(statedefn)
    histogram = zeros(Int, WangLandau.histogram_size(statedefn))
    logdos = zeros(Float64, WangLandau.histogram_size(statedefn))
    logf = 0.1
    new_index = WangLandau.wl_trial!(state, old_index, statedefn, logdos, histogram, logf, ffc)
    @test 1 ≤ new_index ≤ length(histogram)

    # nc = WangLandau.NoCatchup()
    # sim2 = CommonSolve.init(prob; catchup_strategy = nc)
    # @test WangLandau.catchup_enabled(nc) == false
    # @test WangLandau.catchup_value(nc) == 0.0
    # WangLandau.update!(nc, sim2)  
end

@testset "DynamicFractionalCatchup" begin
    Random.seed!(1234)

    L = 5
    statedefn = Ising2D(L; periodic = false)
    prob = WangLandauProblem(statedefn)

    dfc = WangLandau.DynamicFractionalCatchup()
    sim = CommonSolve.init(prob; catchup_strategy = dfc)
    sim.logf_strategy = WangLandau.ReduceByFactor(; final = 1e-6)
    sim.logdos .= 1e-8 .* (1:length(sim.logdos))

    WangLandau.update!(dfc, sim)
    val = WangLandau.catchup_value(dfc)
    @test val isa Float64
    @test val ≥ 0.0
    @test WangLandau.catchup_enabled(dfc) == true

    logdos .= 0.0
    histogram .= 0
    logf = 0.01
    new_index = WangLandau.wl_trial!(state, old_index, statedefn, logdos, histogram, logf, dfc)
    @test 1 ≤ new_index ≤ length(histogram)
end

# @testset "DynamicFractionalCatchup" begin
#     using Random
#     Random.seed!(1234)

#     L = 5
#     statedefn = Ising2D(L; periodic = false)
#     prob = WangLandauProblem(statedefn)

#     sim = CommonSolve.init(prob; check_sweeps = 10,final_logf = 1e-3)

#     @test isa(sim, WangLandau.WangLandauSimulation)
#     @test sim.check_steps == 10 * WangLandau.system_size(statedefn)

#     state, old_index = initialise_state(statedefn)
#     logdos = fill(1e-8, WangLandau.histogram_size(statedefn))
#     histogram = zeros(Int, WangLandau.histogram_size(statedefn))
#     logf = 0.1
#     dfc = WangLandau.DynamicFractionalCatchup()

#     new_index = WangLandau.wl_trial!(state, old_index, statedefn, logdos, histogram, logf, dfc)
#     @test 1 ≤ new_index ≤ length(histogram)
#     @test any(histogram .> 0)
#     @test any(logdos .> 0)

#     hist = zeros(Int, size(sim.samples))
#     task_samples = zeros(Int, max(1, sim.tasks_per_thread * Threads.nthreads()))
#     chunk_size = max(1, sim.check_steps ÷ length(task_samples))

#     flag = CommonSolve.step!(sim, hist, task_samples, chunk_size)
#     @test isa(flag, Bool)
#     @test sim.flat_checks ≥ 1
#     @test sim.total_steps ≥ sim.check_steps

#     sim_short = CommonSolve.init(prob;check_sweeps = 5,final_logf = 1e-2,max_total_steps = 200)
#     sim_done = CommonSolve.solve!(sim_short)
#     @test isa(sim_done, WangLandau.WangLandauSimulation)
#     @test sim_done.total_steps > 0
#     @test sim_done.elapsed_time ≥ 0.0
# end
