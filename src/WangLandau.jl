module WangLandau

using CommonSolve: CommonSolve, init, solve, solve!, step!
using ProgressLogging: @logprogress, @withprogress
using TerminalLoggers: TerminalLogger
using Logging: ConsoleLogger, global_logger
using StaticArrays: StaticArrays, SVector, MVector
using Random: Random, seed!, shuffle!
import TOML

const PACKAGE_VERSION = VersionNumber(TOML.parsefile(pkgdir(@__MODULE__, "Project.toml"))["version"])

@doc """
    WangLandau

This is `WangLandau` $PACKAGE_VERSION.
Read the documentation [online]().
"""
WangLandau

export init, solve!, solve, step!  # from CommonSolve
export WangLandauProblem
export CatchupStrategy, NoCatchup, FixedFractionalCatchup#, DynamicFractionalCatchup
export DosIncrementStrategy, ReduceByFactor
export FlatHistogramStrategy, FractionOfMean, StableNumVisits

export initialise_state, random_trial!, commit_trial!, histogram_size

include("strategies/catchup.jl")
include("strategies/f_increment.jl")
include("strategies/flatness.jl")
include("WangLandauProblem.jl")
include("WangLandauSimulation.jl")

function __init__()
    # setup logging
    if isa(stderr, Base.TTY) && (get(ENV, "CI", nothing) ≠ true) # running in terminal?
        global_logger(TerminalLogger()) # enable progress bar
    else
        global_logger(ConsoleLogger())
    end
    Base.global_logger()
end

end
