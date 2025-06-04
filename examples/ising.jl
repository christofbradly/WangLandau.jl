using WangLandau

# simple 2D Ising model -- no field
struct Ising2D{P,I}
    initial_state::Matrix{Int}
    L::Int
    maxE::Int
end
function Ising2D(L; initial_state = nothing, periodic = false)
    maxE = 2 * (L - !periodic)^2 + 2 * (L - !periodic)  # equiv to: ising_full_energy(ones(Int, L, L))
    if isnothing(initial_state)
        initial_state = zeros(Int, L, L)
        initial_flag = false
    else
        # Retain this feature for debugging
        is_err = ArgumentError("Invalid initial state.")
        initial_state isa Matrix{Int} || throw(is_err)
        size(initial_state) == (L, L) || throw(is_err)
        minimum(initial_state) in (-1,1) || throw(is_err)
        maximum(initial_state) in (-1,1) || throw(is_err)
        0 in initial_state && throw(is_err)
        initial_flag = true
    end
    return Ising2D{periodic,initial_flag}(initial_state, L, maxE)
end

#
# computing energy and histogram index
#
function ising_energy_free(grid)
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

function ising_energy_periodic(grid)
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
function local_energy_free(grid, site)
    L = size(grid)[1]
    i, j = Tuple(site)
    ΔE = 0
    ΔE += i > 1 ? grid[i-1,j] : 0   # above
    ΔE += i < L ? grid[i+1,j] : 0   # below
    ΔE += j > 1 ? grid[i,j-1] : 0   # left
    ΔE += j < L ? grid[i,j+1] : 0   # right
    return ΔE * grid[site]
end

function local_energy_periodic(grid, site)
    L = size(grid)[1]
    i, j = Tuple(site)
    ΔE = 0
    ΔE += i > 1 ? grid[i-1,j] : grid[L,j]   # above
    ΔE += i < L ? grid[i+1,j] : grid[1,j]   # below
    ΔE += j > 1 ? grid[i,j-1] : grid[i,L]   # left
    ΔE += j < L ? grid[i,j+1] : grid[i,1]   # right
    return ΔE * grid[site]
end

#
# WangLandau.jl API
#
WangLandau.histogram_size(statedefn::Ising2D{P}) where {P} = (statedefn.maxE ÷ (1 + P) + 1, )
WangLandau.system_size(statedefn::Ising2D) = statedefn.L^2

function WangLandau.initialise_state(statedefn::Ising2D{Periodic,Init}) where {Periodic,Init}
    L = statedefn.L
    if Init
        spins = copy(statedefn.initial_state)
    else
        spins = rand((-1, 1), (L, L))
    end
    if Periodic
        energy = ising_energy_periodic(spins)
    else
        energy = ising_energy_free(spins)
    end
    Eindex = (energy + statedefn.maxE) ÷ 2 ÷ (1 + Periodic) + 1
    return spins, Eindex
end

function WangLandau.random_trial!(_, statedefn::Ising2D)
    L = statedefn.L
    site = rand(CartesianIndices((L, L)))
    return site
end

function WangLandau.hist_index(spins, statedefn::Ising2D{Periodic}, site, old_index) where {Periodic}
    if Periodic
        ΔE = local_energy_periodic(spins, site) ÷ 2
    else
        ΔE = local_energy_free(spins, site)
    end
    new_index = old_index + ΔE
    energy = Periodic ? ising_energy_periodic(spins) : ising_energy_free(spins)
    actual = (energy + statedefn.maxE) ÷ 2 ÷ (1 + Periodic) + 1
    println((old_index, new_index, actual, energy, spins))
    if !(1 ≤ new_index ≤ statedefn.maxE ÷ (1 + Periodic) + 1)
        println((site, old_index, new_index, ΔE, spins, Periodic, L, statedefn.maxE))
        throw(ErrorException("New energy is invalid."))
    end
    return new_index
end

function WangLandau.commit_trial!(spins, ::Ising2D, site, _, _)
    spins[site] *= -1
    return spins
end
