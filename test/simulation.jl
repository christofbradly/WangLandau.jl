using WangLandau
using Test
using Random
using CommonSolve
using Logging 

include(joinpath(@__DIR__, "..", "examples", "ising.jl"))

@testset "WangLandauSimulation" begin
    Random.seed!(12345)

    L = 5
    periodic = false
    statedefn = Ising2D(L; periodic)

    sim = WangLandau.WangLandauSimulation(statedefn)
    @test isa(sim.logf_strategy, ReduceByFactor)
    @test isa(sim.flat_strategy, FractionOfMean)           
    @test isa(sim.catchup_strategy, NoCatchup)    

    sim_check = WangLandau.WangLandauSimulation(statedefn; check_sweeps = 10)
    @test sim_check.check_steps == 10 * WangLandau.system_size(statedefn)

    sim_steps = WangLandau.WangLandauSimulation(statedefn; max_total_steps = 1)
    @test_logs (:warn, r".*") WangLandau.WangLandauSimulation(statedefn; max_total_steps = 1)

    sim_logf = WangLandau.WangLandauSimulation(statedefn; final_logf = 1e-3)
    @test WangLandau.final_value(sim_logf.logf_strategy) == 1e-3

    fs = FractionOfMean(0.75)
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

    sim = init(prob; check_sweeps = 10, final_logf = 1e-3)

    state, old_index = initialise_state(statedefn)
    logdos = zeros(Float64, WangLandau.histogram_size(statedefn))
    histogram = zeros(Int, WangLandau.histogram_size(statedefn))
    logf = 0.1
    catchup = NoCatchup()

    new_index = WangLandau.wl_trial!(state, old_index, statedefn, logdos, histogram, logf, catchup)
    @test 1 ≤ new_index ≤ length(histogram)
    @test any(histogram .> 0)  
    @test any(logdos .> 0)     

    hist = zeros(Int, size(sim.samples))
    task_samples = zeros(Int, max(1, sim.tasks_per_thread * Threads.nthreads()))
    chunk_size = max(1, sim.check_steps ÷ length(task_samples))

    flag = step!(sim, hist, task_samples, chunk_size)
    @test sim.flat_checks ≥ 1
    @test sim.total_steps ≥ sim.check_steps

    sim_short = init(prob; check_sweeps = 5, final_logf = 1e-2, max_total_steps = 200)
    sim_done = solve!(sim_short)
    @test isa(sim_done, WangLandau.WangLandauSimulation)
    @test sim_done.total_steps > 0
    @test sim_done.elapsed_time ≥ 0.0

    state2, old_index2 = initialise_state(statedefn)
    hist2 = zeros(Int, WangLandau.histogram_size(statedefn))
    logdos2 = zeros(Float64, WangLandau.histogram_size(statedefn))
    logf2 = 0.1
    ffc = FixedFractionalCatchup(0.25)

    new_index2 = WangLandau.wl_trial!(state2, old_index2, statedefn, logdos2, hist2, logf2, ffc)
    @test 1 ≤ new_index2 ≤ length(hist2)
end

@testset "FixedFractionalCatchup" begin
    Random.seed!(12345)

    L = 5
    periodic = false
    statedefn = Ising2D(L; periodic)

    ffc = FixedFractionalCatchup(0.25)
    sim_stub = (; logdos = [1e-8, 2e-8, 3e-8, 4e-8], catchup_strategy = ffc)

    oldmin = ffc.minval
    WangLandau.update!(ffc, sim_stub)
    @test ffc.minval == oldmin
    @test ffc.minval == 1e-8  
    @test WangLandau.catchup_enabled(ffc) == true
    @test WangLandau.catchup_value(ffc) ≈ ffc.minval * ffc.fraction

    ffc_zero = FixedFractionalCatchup(0.25)
    sim_zero = (; logdos = zeros(4), catchup_strategy = ffc_zero)
    WangLandau.update!(ffc_zero, sim_zero)
    @test ffc_zero.minval == 0.0       
end

@testset "FlatHistogramStrategy" begin
    Random.seed!(1234)

    @testset "FractionOfMean" begin
        f1 = FractionOfMean(0.8)

        @test_throws ArgumentError FractionOfMean(1.2)
        @test_throws ArgumentError FractionOfMean(-0.5)
        @test_throws ArgumentError FractionOfMean(0.5, -1)

        f = FractionOfMean(0.8)
        hist = zeros(Int, 5)
        @test_throws ErrorException WangLandau.isflat(f, hist)

        hist .= [0, 0, 3, 0, 0]
        @test !WangLandau.isflat(f, hist)

        hist .= [5, 5, 5, 5, 5]
        @test WangLandau.isflat(f, hist)

        hist .= [10, 9, 7, 6, 1]
        @test !WangLandau.isflat(f, hist)

        flat_strategy = FractionOfMean(0.8)
        sim_stub = (; flat_iterations = 1, total_steps = 10)
        @test isnothing(WangLandau.update!(flat_strategy, sim_stub))
        @test isnothing(WangLandau.update!(flat_strategy, (;)))   
    end

    @testset "StableNumVisits" begin
        @test_throws ArgumentError StableNumVisits(2, 0)

        s1 = StableNumVisits(1, 10)
        hist = [0, 1, 2, 3, 4]
        @test !WangLandau.isflat(s1, hist)
        @test s1.stablesteps == 0

        s2 = StableNumVisits(2, 5)
        hist .= [2, 2, 2, 2, 2]
        @test !WangLandau.isflat(s2, hist)    
        @test s2.numvisits == 0               
        @test s2.stablesteps == 0             
        s2.stablesteps = 6                    
        @test WangLandau.isflat(s2, hist)

        sim_stub = (; flat_iterations = 2, total_steps = 50)
        s3 = StableNumVisits(2, 5)
        WangLandau.update!(s3, sim_stub)
        @test s3.numvisits == 0
        @test s3.stablesteps == 0
    end
end