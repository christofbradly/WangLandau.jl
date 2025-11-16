using WangLandau
using BenchmarkTools
using Random

include(joinpath(@__DIR__, "..", "examples", "ising.jl"))

function make_sim(L::Int; kwargs...)
    Random.seed!(12345)

    periodic = false
    statedefn = Ising2D(L; periodic)
    problem = WangLandauProblem(statedefn)

    sim = solve(problem;
        check_sweeps=check_sweeps,
        final_logf=final_logf,
        tasks_per_thread=tasks_per_thread,
    )

    return sim
end

@info "Benchmark running with $(Threads.nthreads()) threads"

const SUITE = @benchmarkset "WangLandau" begin
    @benchmarkset "Scaling" begin
        @case "L=16" begin
            make_sim(16; check_sweeps=500, final_logf=1e-3)
        end seconds=60

        @case "L=32" begin
            make_sim(32; check_sweeps=1000, final_logf=1e-3)
        end seconds=60
    end

    @benchmarkset "Threading" begin
        @case "L=16" begin
            make_sim(16; check_sweeps=1000, final_logf=1e-3)
        end seconds=60

        @case "L=32" begin
            make_sim(32; check_sweeps=1000, final_logf=1e-3)
        end seconds=60
    end
end