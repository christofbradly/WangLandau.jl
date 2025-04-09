# Advanced Usage

## Strategies

The Wang-Landau algorithm has rules for determining the histogram flatness and how to update ``f``. These can be controlled by via keyword arguments ``flat_strategy`` and ``logf_strategy`` passed to `solve`.

### Histogram flatness criterion

The flatness criterion is defined by [`FlatHistogramStrategy`](@ref). The default `FlatHistogramStrategy` is `FractionOfMean` which determines flatness by whether the minimum value of the histogram is at least a fraction `x` of the mean value. Because it is the default strategy, the parameter `x` can also be supplied directly to `solve(prob; kwargs...)` with the keyword argument `flat_tolerance = x`.
A custom `strat::FlatHistogramStrategy` can be passed to the algorithm via `solve(prob; flat_strategy = strat)`.

Implemented flatness strategies include
- [`FractionOfMean`](@ref): 
- [`StableNumVisits`](@ref): 

### Density of states parameter ``f``

The alteration of ``f`` (actually ``\\log f``) at each iteration is defined by [`DosIncrementStrategy`](@ref).
The default `DosIncrementStrategy` is `ReduceByFactor` which defines a parameter `x` (default `x=0.5`) and updates ``\\log f \\to x\\log f`` after each flatness iteration. It also defines initial and final values for ``\\log f`` with defaults `1.0` and `1e-6`, respectively. Because it is the default strategy, the final value of ``\\log f`` can also be supplied directly to `solve(prob; kwargs...)` with the keyword argument `final_logf`.
A custom `strat::DosIncrementStrategy` can be passed to the algorithm via `solve(prob; logf_strategy = strat)`.

Implemented ``\\log f`` strategies include
- [`ReduceByFactor`](@ref):

## Other features

- parallelisation (multi-threading): implemented with `@atomic`. Prior to Julia v.1.12 this requires the [`Atomix.jl`](https://github.com/JuliaConcurrent/Atomix.jl) package

- parallelisation (replica exchange): TBC
- multi-index histogram/dos: Typically the energy ``E`` is a well-defined physical property of the system being studied, but internally `WangLandau.jl` only treats it as an index to the histogram and density of states. Thus, one could have a system with two parameters of interest, e.g. energy and order parameter. See [`random_trial!`](@ref)
- user defined `DefType` and `StateType`
