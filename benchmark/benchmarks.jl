using WangLandau
using BenchmarkTools

include(joinpath(@__DIR__, "..", "examples", "ising.jl"))

function make_sim(L::Int; kwargs...)
    periodic = false
    statedefn = Ising2D(L; periodic)
    problem = WangLandauProblem(statedefn)
    return solve(problem; kwargs...)
end

const SUITE = @benchmarkset "WangLandau" begin
    @benchmarkset "Scaling" begin
        @case "L=8"  make_sim(8;  check_sweeps=10_000, final_logf=1e-6)  seconds=10
        @case "L=16" make_sim(16; check_sweeps=10_000, final_logf=1e-6) seconds=10
        @case "L=24" make_sim(24; check_sweeps=10_000, final_logf=1e-6) seconds=10
        @case "L=32" make_sim(32; check_sweeps=10_000, final_logf=1e-6) seconds=10
    end

    @benchmarkset "Threading" begin
        @case "L=16, threads=$(Threads.nthreads())" make_sim(16; check_sweeps=10_000, final_logf=1e-6) seconds=15
    end
end