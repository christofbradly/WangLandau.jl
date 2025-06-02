using WangLandau

# simple 2D Ising model - free boundary conditions, no field

function ising_free_energy(grid)
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

function ising_periodic_energy(grid)
    L = size(grid)[1]
    energy = 0
    for i in 1:L
        for j in 1:L
            energy += grid[i,j] * grid[mod1(i + 1, L),j]   # below
            energy += grid[i,j] * grid[i,mod1(j + 1, L)]   # right
        end
    end

    return -energy
end

# actual energy change is twice local energy
function local_free_energy(grid, site)
    L = size(grid)[1]
    i, j = Tuple(site)
    ΔE = 0
    ΔE += i > 1 ? grid[i-1,j] : 0   # above
    ΔE += i < L ? grid[i+1,j] : 0   # below
    ΔE += j > 1 ? grid[i,j-1] : 0   # left
    ΔE += j < L ? grid[i,j+1] : 0   # right
    return ΔE * grid[site]
end

function local_periodic_energy(grid, site)
    L = size(grid)[1]
    i, j = Tuple(site)
    ΔE = 0
    ΔE += i > 1 ? grid[i-1,j] : grid[L,j]   # above
    ΔE += i < L ? grid[i+1,j] : grid[1,j]   # below
    ΔE += j > 1 ? grid[i,j-1] : grid[i,L]   # left
    ΔE += j < L ? grid[i,j+1] : grid[i,1]   # right
    return ΔE * grid[site]
end

mutable struct Ising2D{P}
    spins::Matrix{Int}
    const L::Int
    const maxE::Int
    energy::Int
end
function Ising2D(L; initial_state = nothing, periodic = false)
    maxE = 2 * (L - !periodic)^2 + 2 * (L - !periodic)  # equiv to: ising_full_energy(ones(Int, L, L))
    if isnothing(initial_state)
        initial_state = rand((-1, 1), (L, L))
    else
        is_err = ArgumentError("Invalid initial state.")
        initial_state isa Matrix{Int} || throw(is_err)
        size(initial_state) == (L, L) || throw(is_err)
        minimum(initial_state) in (-1,1) || throw(is_err)
        maximum(initial_state) in (-1,1) || throw(is_err)
        0 in initial_state && throw(is_err)
    end
    if periodic
        E0 = ising_periodic_energy(initial_state)
    else
        E0 = ising_free_energy(initial_state)
    end
    Eindex = (E0 + maxE) ÷ 2 ÷ (1 + periodic) + 1
    return Ising2D{periodic}(L, maxE, Eindex, initial_state)
end

function Base.copy(state::Ising2D{P}) where {P}
    return Ising2D{P}(state.L, state.maxE, state.energy, copy(state.spins))
end

# WangLandau.jl API
WangLandau.histogram_size(state::Ising2D{P}) where {P} = (state.maxE ÷ (1 + P) + 1, )
WangLandau.system_size(state::Ising2D) = state.L^2

function WangLandau.initialise_state(state::Ising2D{P}) where {P}
    return Ising2D(state.L; periodic = P)
end

function WangLandau.random_trial!(state::Ising2D{Periodic}) where {Periodic}
    L = state.L
    site = rand(CartesianIndices((L, L)))
    oldE = state.energy
    ΔE = Periodic ? local_periodic_energy(state.spins, site) ÷ 2 : local_free_energy(state.spins, site)
    newE = oldE + ΔE
    return site, oldE, newE
end

function WangLandau.commit_trial!(state::Ising2D, site, _, newE)
    state.spins[site] *= -1
    state.energy = newE
    1 ≤ state.energy ≤ 2state.maxE + 1 || throw(ErrorException("New energy is invalid."))
    return state
end
