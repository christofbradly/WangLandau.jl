#!/usr/bin/env julia
using JSON, Printf

root = ARGS[1]
latest = get(ENV, "LATEST_TAG", "UNKNOWN")

println("## Benchmark Results\n")

# Look for both thread configurations
for threads in ["1", "4"]
    dir = joinpath(root, "results-julia-1.11-threads-$threads", "julia-1.11", "threads-$threads")

    basefile = joinpath(dir, "results_WangLandau@$latest.json")
    headfile = joinpath(dir, "results_WangLandau@HEAD.json")

    if !isfile(basefile) || !isfile(headfile)
        println("⚠ Missing files for threads=$threads — skipping.\n")
        continue
    end

    baseline = JSON.parsefile(basefile)
    current  = JSON.parsefile(headfile)

    println("<details><summary>Threads: $threads</summary>\n")
    println("| Benchmark | $latest | HEAD | Ratio |")
    println("|-----------|--------|------|-------|")

    for key in sort(collect(keys(baseline)))
        old = baseline[key]
        new = current[key]
        ratio = new / old
        @printf("| %s | %.5g | %.5g | %.3f |\n", key, old, new, ratio)
    end

    println("\n</details>\n")
end