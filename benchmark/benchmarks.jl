using WangLandau
using CommonSolve
using BenchmarkTools

include(joinpath(@__DIR__, "..", "examples", "ising.jl"))

function make_sim(L::Int; kwargs...)
    statedefn = WangLandau.Ising2D(L)
    problem = WangLandau.WangLandauProblem(statedefn)
    sim = CommonSolve.init(problem; kwargs...)
    return sim
end

const SUITE = @benchmarkset "WangLandau" begin
    @benchmarkset "Small " begin
        @case "L=8" begin
            sim = make_sim(8; check_steps=5_000, max_total_steps=50_000)
            CommonSolve.solve!(sim)
        end seconds=10

        @case "L=16" begin
            sim = make_sim(16; check_steps=10_000, max_total_steps=100_000)
            CommonSolve.solve!(sim)
        end seconds=15
    end

    @benchmarkset "Scaling" begin
        for L in (8, 16, 24, 32)
            @case "L=$L" begin
                local L_ = $L
                sim = make_sim(L; check_steps=10_000, max_total_steps=100_000)
                CommonSolve.solve!(sim)
            end seconds=10
        end
    end

    @benchmarkset "Threading" begin
        @case "L=16, threads=$(Threads.nthreads())" begin
            sim = make_sim(16; check_steps=10_000, max_total_steps=100_000)
            CommonSolve.solve!(sim)
        end seconds=15
    end
end