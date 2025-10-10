using Test
using WangLandau

@testset "WangLandauProblem" begin
    include("problem.jl")
end

@testset "WangLandauSimulation" begin
    include("simulation.jl")
end