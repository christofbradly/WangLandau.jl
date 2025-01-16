module WangLandau

using CommonSolve: CommonSolve, init, solve, solve!, step!
using ProgressLogging: ProgressLogging, @logprogress, @withprogress
using TerminalLoggers: TerminalLoggers
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
public initialise_state, random_move, test_move, commit!, histogram_size
export NoCatchup, FixedFractionalCatchup
export LogReduceByFactor

include("WangLandauProblem.jl")
include("WangLandauSimulation.jl")
include("strategies/catchup.jl")
include("strategies/f_increment.jl")
include("strategies/flatness.jl")

function __init__()
    # setup logging
end

end
