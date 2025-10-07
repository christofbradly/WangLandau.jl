using WangLandau
using Test
using Random

include(joinpath(@__DIR__, "..", "examples", "ising.jl"))

@testset "WangLandauProblem" begin
    Random.seed!(12345)

    # Define Problem
    L = 5
    periodic = false
    statedefn = Ising2D(L; periodic)
    prob = WangLandauProblem(statedefn)

    @test WangLandau.histogram_size(statedefn) == (41,)
    @test WangLandau.system_size(statedefn) == 25

    state, index = WangLandau.initialise_state(statedefn)
    @test state == [
         1  -1   1  -1  -1;
         1   1  -1  -1   1;
        -1  -1  -1   1  -1;
        -1  -1  -1   1   1;
        -1   1   1  -1  -1
    ]
    @test index == 23

    trial = random_trial!(state, prob.statedefn)
    @test trial == CartesianIndex(4, 1)

    new_index = histogram_index(state, statedefn, trial, index)
    @test new_index == 26

    commmit_state = WangLandau.commit_trial!(copy(state), statedefn, trial, index, new_index)
    @test commmit_state == [
        1  -1   1  -1  -1;
        1   1  -1  -1   1;
       -1  -1  -1   1  -1;
        1  -1  -1   1   1;
       -1   1   1  -1  -1
    ]

    revert_state = WangLandau.revert_trial!(copy(commmit_state), statedefn, trial, new_index, index)
    @test revert_state == commmit_state
end