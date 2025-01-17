using WangLandau

# simple 2D Ising model - free boundary conditions, no field

function ising_full_energy(grid)
    L = size(grid)[1]
    energy = 0
    for i in 1:L-1
        for j in 1:L-1
            energy += grid[i,j] * grid[i+1,j]   # below
            energy += grid[i,j] * grid[i,j+1]   # right
        end
        energy += grid[L,i] * grid[L,i+1]   # last row
        energy += grid[i,L] * grid[i+1,L]   # last column
    end

    return -energy
end

function local_energy(grid, L, site)
    i, j = Tuple(site)
    ΔE = 0
    ΔE += i > 1 ? grid[i-1,j] : 0   # above
    ΔE += i < L ? grid[i+1,j] : 0   # below
    ΔE += j > 1 ? grid[i,j-1] : 0   # left
    ΔE += j < L ? grid[i,j+1] : 0   # right
    return ΔE
end

mutable struct Ising2D
    L::Int
    maxE::Int
    energy::Int
    spins::Matrix{Int}
end
function Ising2D(L)
    maxE = 2 * (L - 1)^2 + 2 * (L - 1)  # equiv to: ising_full_energy(ones(Int,L,L))
    spins = rand((-1, 1), (L, L))
    E0 = ising_full_energy(spins) + maxE + 1
    return Ising2D(L, maxE, E0, spins)
end

# WangLandau.jl API
WangLandau.initialise_state(state::Ising2D) = state

WangLandau.histogram_size(state::Ising2D) = (2state.maxE + 1, )

function WangLandau.random_move(state::Ising2D)
    L = state.L
    return rand(CartesianIndices((L, L)))
end

function WangLandau.test_move(state::Ising2D, site)
    ΔE = state.spins[site] * local_energy(state.spins, state.L, site)
    return state.energy + ΔE
end

function WangLandau.commit!(state::Ising2D, site, newE)
    state.spins[site] *= -1
    state.energy = newE
    1 ≤ state.energy ≤ 2state.maxE + 1 || throw(ErrorException("New energy is invalid."))
    return state
end

# Run simulation
L = 10
prob = WangLandauProblem(Ising2D(L))

sim = solve(prob)
