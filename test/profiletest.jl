include("../examples/ising.jl")

L = 10
prob = WangLandauProblem(Ising2D(L; periodic = true))
solve(prob; final_logf = 1e-1)
@profview solve(prob)
