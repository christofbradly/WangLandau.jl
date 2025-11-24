#!/usr/bin/env julia

using JSON
using AirspeedVelocity

root = ARGS[1]
latest = get(ENV, "LATEST_TAG", "unknown")

println("## Benchmark Results\n")

for dir in sort(readdir(root))
    full = joinpath(root, dir)

    # Expect directories like: results-julia-1.11-threads-1
    startswith(dir, "results-julia") || continue

    # Extract thread count
    threads = try
        parse(Int, split(dir, "-")[end])
    catch
        println("⚠ Could not extract thread count from $dir — skipping.")
        continue
    end

    # Real json path (note the nested “results” folder!)
    jsondir = joinpath(full, "results", "julia-1.11", "threads-$threads")

    if !isdir(jsondir)
        println("⚠ Missing files for threads=$threads — skipping.\n")
        continue
    end

    println("<details><summary>Threads: $threads</summary>\n")

    try
        table = AirspeedVelocity.benchpkgtable(
            "",
            input_dir=jsondir,
            rev="$latest,HEAD",
            mode="ratio"
        )
        println(table)
    catch e
        println("❌ Error generating table for threads=$threads: ", e)
    end

    println("\n</details>\n")
end
