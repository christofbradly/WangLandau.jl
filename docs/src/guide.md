# WangLandau.jl

WangLandau.jl is a pure Julia implementation of the Wang-Landau Monte Carlo algorithm for statistical mechanics problems. The package defines a flexible interface that is designed to apply the algorithm to any user-defined problem.

## Wang-Landau algorithm

The Wang-Landau algorithm works by randomly sampling from a state space and iteratively updating a corresponding density of states until some accuracy criterion is achieved. 

Suppose we initialise a system in state ``i`` (e.g. an Ising configuration on a lattice), which has energy ``E_i``, and we propose a trial move to state ``j``, which has energy ``E_j``.
The proposed state is accepted based on a Metropolis-like condition
```
P(i \\to j) = \\min(1, \\frac{g(E_i)}{g(E_j)})
```
where ``g(E)`` is the density of states at energy ``E``. Classically, the density of states is motivated by the Boltzmann distribution ``g(E) = \\exp(-E)``. A random roll is performed against ``P(i \\to j)`` and either the proposed state or the original state is accepted as the next state.
The sampled states are counted with a histogram ``h(E)``, initially set to 0.
When the next state, with energy ``E``, is accepted, ``h`` is incremented to count the new state: ``h(E) \\to h(E) + 1`` and the density of states is increased according to ``g(E) \\to  f * g(E)``, where ``f`` is a dynamic parameter.

A single iteration of the Wang-Landau algorithm is as follows: After a set number of trial moves, the histogram is checked to see if it is sufficiently 'flat', for example, ``\\min h(E) \\geq 0.8\\bar{h(E)}``, where ``\\bar{h(E)}`` is the average value of ``h``. If the histogram is flat, then ``h`` is reset to 0, the parameter ``f`` is decreased, for example ``f \\to \\sqrt{f}``, and a new iteration commences. 

The algorithm terminates when ``f`` reaches a predetermined size, for example, one may initialise ``f = e\approx 2.718\\ldots`` and after 14 iterations ``f < 10^{-6}``. At this point the error in the density of states is ``\sqrt{f}``, and thermodynamic properties can be calculated from ``g(E)``. For numerical accuracy, the density of states and the parameter ``f`` are usually handled with logarithms.

## Interface

WangLandau.jl adheres to the ['CommonSolve'](https://github.com/SciML/CommonSolve.jl) interface.
```julia
prob = WangLandauProblem(::StateType)
sim = CommonSolve.solve(prob; kwargs...)
```
The `StateType` should be defined by the user for a specific problem. Several methods should be defined by the user to dispatch on `StateType`, see [`WangLandauProblem`](@ref) for details.
The `kwargs` passed to `solve` define the parameters of the algorithm, and are ultimately passed to [`WangLandauSimulation`](@ref), which is the main internal `struct` and the output type of `sim`. `WangLandauSimulation` is not exported and should not be constructed explicitly, it is handled via the `CommonSolve` interface.

The main parameters include:
- `check_sweeps`: A statistical mechanics problem will have some finite size `N` and a Monte Carlo sweep is `N` trial moves. `check_sweeps` controls how many sweeps to perform before checking if the histogram is flat. Checking too often can be inefficient and checking too rarely will mean the algorithm runs longer than necessary.
- `flat_tolerance`: defines the number `x` such that ``h(E)`` is considered flat when ``\\min h \\geq x \\bar{h}``
- `final_logf`: when ``f`` reaches this value the algorithm terminates

The output `sim::WangLandauSimulation` will contain two `Array`s:
- `sim.samples`: the histogram of samples ``h(E)``
- `sim.logdos`: The (logarithm) of the density of states ``\\log g(E)``

## References
- Original paper: [Wang & Landau, *PRL* 2001](https://link.aps.org/doi/10.1103/PhysRevLett.86.2050)
