#
# update path to local copy of repo
#
include("C:/University/Research/WangLandau.jl/examples/ising.jl")

# define problem
L = 10
periodic = false
statedefn = Ising2D(L; periodic)

# Run simulation
prob = WangLandauProblem(statedefn)
sim = solve(prob; check_sweeps = 100, final_logf = 1e-6)

##########################################################
using Plots

# extract data from simulation
Erange = -(sim.statedefn.maxE):(2*(1 + periodic)):(sim.statedefn.maxE)
mask = sim.samples .> 0     # check for unsampled states
avg = sum(sim.samples[mask]) / length(sim.samples[mask])
std = sqrt(sum((sim.samples[mask] .- avg) .^ 2) / (length(mask) - 1))
logdos_shift = maximum(sim.logdos[mask])
dos = exp.(sim.logdos[mask] .- logdos_shift)

# Plot histogram (should be flat!)
plt1 = scatter(Erange, sim.samples; label = :none, xlabel = "E", ylabel = "samples", msw = 0, ms = 2)
hline!(plt1, [avg, sim.flat_strategy.tol * avg]; c = :grey, ls = :dash, label = :none)
scatter!(plt1, [0], [avg], err=[std]; label = :none, ms=1,msw=1, c = :black)

# plot density of states
plt2 = scatter(Erange[mask], dos; label = :none, xlabel = "E", ylabel = "g(E)", yscale = :log10, ms = 2, msw = 0)

plot(plt1, plt2; layout=2)
