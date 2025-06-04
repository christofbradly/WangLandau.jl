```@meta
ShareDefaultModule=true
```
# Example: 2D Ising model

```@contents
Pages = ["ising.md"]
Depth = 2:3
```

## The Ising model
The Ising model is a classic statistical mechanics model defined by a Hamiltonian
```math
H = -J \sum_{i \sim j} S_i S_j
```
where ``S_i \in (-1,1)`` is the spin of site ``i``. The sum is over all adjacent pairs on a specified lattice. Here we assume that the interaction strength is in units where ``J = 1`` and we ignore any external magnetic field. We will apply this model on a square grid of size ``L \times L`` sites. This model is the simplest statistical mechanical system that exhibits a continuous phase transition.


## Setup
To apply the Wang-Landau algorithm, we load `WangLandau` along with necessary definitions for the Ising model
```@example
using WangLandau, LinearAlgebra, Plots, LaTeXStrings;

include("../../examples/ising.jl");
```

The main definition in `ising.jl` is the mutable struct `Ising2D`, which holds an integer array representing the configuration of spins and the energy of the configuration. The rest of the code defines methods for applying the `WangLandau` interface to `Ising2D` (see [`WangLandauProblem`](@ref)). These methods define the histograms and specify how to make and accept random trial moves. For the Ising model a simple trial move is flipping a single spin at a random site on the lattice.

This example can have periodic or free boundary conditions, here we will define the problem for the case of free boundary conditions:
```@example
L = 15;

P = false;

prob = WangLandauProblem(Ising2D(L; periodic=P));
```

When `Ising2D` is constructed it will initialise a random configuration of spins.

## Simulation

Next, we run the simulation. We will use the default strategies for altering the ``f`` parameter, but set a final desired value of ``f < 10^{-6}``. We can also set how often to check for flatness by setting the keyword `check_sweeps`. A single Monte Carlo sweep is ``N`` trial moves, where ``N`` is the canonical size of the system, in this case ``N = L^2``.
```@example
sim = solve(prob; check_sweeps = 500, final_logf = 1e-6);
```

To see if the algorithm performed well, we should inspect the histogram of samples. The output `sim.samples` is the accumulation of the histogram at each iteration. 
```@example
scatter(sim.samples; label = :none, ylabel = "samples", msw = 0, ms = 2)
```
If the histogram is flat then we can proceed to analysis.

The samples histogram is indexed by ``i \in 1,\ldots,i_\mathrm{max}``, where ``i_\mathrm{max} = 2 (L - P)^2 + 2 (L - P)`` is the number of accessible energy values and ``P=0`` for free boundary and ``P=1`` for periodic boundary. This is more efficient for internal processing and the conversion to consecutive array indexing is handled by the `Ising2D` code. However, for physical analysis we want a range of real energy values. The extremal energy values for an ``L \times L`` lattice are ``\pm E_\mathrm{max}`` where ``E_\mathrm{max} = 2(L - (1 - P))^2 + 2(L - (1 - P))``, and the change in energy from flipping a single spin is ``2(1+P)``.
In fact, ``E_\mathrm{max}`` is already stored in `Ising2D` so we can just define
```@example
maxE = sim.statedefn.maxE;

Es = -maxE:(2*(1 + P)):maxE;
```

Next, we look at the density of states, normalised to its maximum value. We also apply a mask to avoid any hard-to-reach states that were not sampled.
```@example
mask = sim.samples .> 0;

logdos_shift = maximum(sim.logdos[mask]);

dos = exp.(sim.logdos[mask] .- logdos_shift);

scatter(Es[mask], dos; label = :none, xlabel = L"E", ylabel = L"g(E)", yscale = :log10, ms = 2, msw = 0)
```

## Critical point

The 2D Ising model has a continuous phase transition at ``\beta_c J = \log (1+\sqrt{2}) / 2 \approx 0.44\ldots``, where ``\beta`` is inverse temperature. The output of the simulation is the entropy ``\log g(E)``, i.e. the logarithm of the density of states, from which we can calculate moments of the energy
```math
\langle E^n \rangle = \frac{\sum_E E^n g(E) \exp(-\beta E)}{\sum_E g(E) \exp(-\beta E)}
```
To compute moments of the density of states from the output data we define 
```@example
function moment_E(Es, g, n, β); E_factor = @. Es^n * exp(-Es * β); return dot(g, E_factor); end;
```
where `Es` is the range of energy values, `g` is the corresponding density of states, `n` is the desired moment and `β` is an inverse temperature value.
Using this we compute the partition function, the average energy and the specific heat (energy variance) as the 0th, 1st and 2nd moments, for a range of values of ``\beta``
```@example
βs = 0.1:0.01:1.0;

Z = [moment_E(Es[mask], dos, 0, β) for β in βs];

avgE = [moment_E(Es[mask], dos, 1, β) for β in βs];

avgE2 = [moment_E(Es[mask], dos, 2, β) for β in βs];

varE = @. (avgE2 / Z - (avgE / Z)^2) * βs^2;
```

Now we can plot the average energy
```@example
plot(βs, avgE ./ Z ./ maxE; label = :none, xlabel = L"\beta", ylabel = L"\langle E \rangle / E_\mathrm{max}", xlim = (0, Inf))
```
and the specific heat, or variance
```@example
plot(βs, varE ./ L^2; label = :none, xlabel = L"\beta", ylabel = L"\mathrm{var}(E)", xlim = (0, Inf), ylim = (0, Inf))
```

Finally, the peak of the specific heat is a signature of the location of the critical point
```@example
varE_max, varE_max_i = findmax(varE);

varE_max / L^2, βs[varE_max_i]
```
which compares well to the known value.
