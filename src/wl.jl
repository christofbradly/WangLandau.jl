# methods for running WangLandau
L = 10

state = initialise(Ising2D(L))
dims = histogram_size(state)

logdos = zeros(dims)
samples = zeros(Int, dims)
temp_hist = zeros(Int, dims)

logf = 1.0
final_logf = 1e-6
logf_decrease = 0.5
flatness = 0.9
enable_catchup = false
catchup_fraction = 0.95
flat_iterations = 0

check_steps = 100   # number of steps between flatness checks
max_total_steps = 1e6   # global max steps to prevent runaway

while logf > finalf
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
        temp_hist[newE] += 1
        if iszero(logdos[newE])
            if enable_catchup
                logdos[newE] = minimum(logdos[logdos .> 0]) * catchup_fraction
            else
                logdos[newE] += logf
            end
        end
    end
    
    if isflat(temp_hist, flatness)
        logf *= logf_decrease
        flat_iterations += 1
        samples .+= temp_hist
        temp_hist = zeros(Int, dims)
    end

    total_steps += check_steps
    if total_steps > max_total_steps 
        @info "Global maximum number of steps reached, ending simulation."
        break
    end
end

function isflat(samples, x)
    nonzeros = count(samples > 0)
    # avoid false positive during initial warmup
    nonzeros < 0.5 * length(samples) && return false
    return sum(samples) ./ nonzeros > x
end
