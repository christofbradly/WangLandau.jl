using WangLandau
using BenchmarkTools

include(joinpath(@__DIR__, "..", "examples", "ising.jl"))

function make_sim(L::Int; kwargs...)
    periodic = false
    statedefn = Ising2D(L; periodic)
    problem = WangLandauProblem(statedefn)
    return solve(problem; kwargs...)
end

@info "Benchmark running with $(Threads.nthreads()) threads"

const SUITE = @benchmarkset "WangLandau" begin
    @benchmarkset "Scaling" begin
       @case "L=8" begin
            make_sim(32; check_sweeps=500, final_logf=1e-3)
        end seconds=60
        @case "L=16" begin
            make_sim(64; check_sweeps=1_000, final_logf=1e-3)
        end seconds=60
    end

    @benchmarkset "Threading" begin
        @case "L=8" begin
            make_sim(32; check_sweeps=1_000, final_logf=1e-3)
        end seconds=60

        @case "L=16" begin
            make_sim(64; check_sweeps=1_000, final_logf=1e-3)
        end seconds=60
    end
end