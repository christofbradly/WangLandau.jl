module WangLandau

using CommonSolve: CommonSolve, init, solve, solve!, step!
using ProgressLogging: @logprogress, @withprogress
using TerminalLoggers: TerminalLogger
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
export DosIncrementStrategy, LogReduceByFactor
export FlatHistogramStrategy, FractionOfMean, StableNumVisits

public initialise_state, random_move, test_move, commit!, histogram_size, measure

include("strategies/catchup.jl")
include("strategies/f_increment.jl")
include("strategies/flatness.jl")
include("WangLandauProblem.jl")
include("WangLandauSimulation.jl")

function __init__()
    # setup logging
    Base.global_logger(TerminalLogger(right_justify=120))
end

end
