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
       @case "L=8" begin
            make_sim(8; check_sweeps=500, max_total_steps=5_000, final_logf=1e-3)
        end seconds=10
        @case "L=16" begin
            make_sim(16; check_sweeps=1_000, max_total_steps=10_000, final_logf=1e-3)
        end seconds=10
    end

    @benchmarkset "Threading" begin
        @case "L=8, threads=$(Threads.nthreads())" begin
            make_sim(8; check_sweeps=1_000, max_total_steps=10_000, final_logf=1e-3)
        end seconds=10

        @case "L=16, threads=$(Threads.nthreads())" begin
            make_sim(16; check_sweeps=1_000, max_total_steps=10_000, final_logf=1e-3)
        end seconds=10
    end
end