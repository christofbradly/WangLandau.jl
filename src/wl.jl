# methods for running WangLandau

function isflat(hist, tol)
    numnonzeros = length(hist[hist .> 0])
    # avoid false positive during initial warmup
    numnonzeros == 0 && throw(ErrorException("No samples in local histogram."))
    numnonzeros == 1 && return false
    mean = sum(hist[hist .> 0]) / numnonzeros
    gap = minimum(hist[hist .> 0]) / mean
    return gap > tol
end

include("ising.jl")
L = 10

state = initialise(Ising2D(L))
dims = histogram_size(state)

logdos = zeros(dims)
samples = zeros(Int, dims)
global temp_hist = zeros(Int, dims)

global logf = 1.0
final_logf = 1e-6
logf_decrease = 0.5
flatness = 0.9
enable_catchup = false
catchup_fraction = 0.95

check_steps = 100   # number of steps between flatness checks
max_total_steps = 1e8   # global max steps to prevent runaway

# Diagnostics
flat_checks = 0
flat_iterations = 0
total_steps = 0

while logf > final_logf
    for _ in 1:check_steps
        oldE = state.energy
        move = random_move(state)
        newE = test_move(state, move)
        
        olddos = logdos[oldE]
        newdos = logdos[newE]

        if log(rand()) < olddos - newdos
            commit!(state, move, newE)
        else
            newE = oldE
        end
        global temp_hist[newE] += 1
        if iszero(logdos[newE])
            if enable_catchup && flat_iterations > 0
                global logdos[newE] = minimum(logdos[logdos .> 0]) * catchup_fraction
            else
                global logdos[newE] += logf
            end
        else
            global logdos[newE] += logf
        end
    end
    
    if isflat(temp_hist, flatness)
        global logf *= logf_decrease
        global flat_iterations += 1
        samples .+= temp_hist
        global temp_hist = zeros(Int, dims)
    end

    global flat_checks += 1
    global total_steps += check_steps
    if total_steps > max_total_steps 
        @info "Global maximum number of steps reached, ending simulation."
        break
    end
end
