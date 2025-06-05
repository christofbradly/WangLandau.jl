module WangLandau

using CommonSolve: CommonSolve, init, solve, solve!, step!
using Atomix: Atomix, @atomic
using ProgressLogging: @logprogress, @withprogress
using TerminalLoggers: TerminalLogger
using Logging: ConsoleLogger, global_logger
using StaticArrays: StaticArrays, SVector, MVector
using Random: Random, seed!, shuffle!
using Printf: @printf
import TOML

const PACKAGE_VERSION = VersionNumber(TOML.parsefile(pkgdir(@__MODULE__, "Project.toml"))["version"])

@doc """
    WangLandau

This is `WangLandau` $PACKAGE_VERSION.
"""
WangLandau

export init, solve!, solve, step!  # from CommonSolve
export WangLandauProblem
public WangLandauSimulation
export CatchupStrategy, NoCatchup, FixedFractionalCatchup#, DynamicFractionalCatchup
export DosIncrementStrategy, ReduceByFactor
export FlatHistogramStrategy, FractionOfMean, StableNumVisits

export initialise_state, histogram_size, system_size
export random_trial!, histogram_index, commit_trial!, revert_trial!

include("strategies/f_increment.jl")
include("strategies/flatness.jl")
include("strategies/catchup.jl")
include("WangLandauProblem.jl")
include("WangLandauSimulation.jl")

function __init__()
    # setup logging - from Rimu.jl
    if isdefined(Main, :IJulia) && Main.IJulia.inited # are we running in Jupyter?
        # need for now as TerminalLoggers currently does not play nice with Jupyter
        # install a bridge to use ProgressMeter under the hood
        # may become unneccessary in the future
        ConsoleProgressMonitor.install_logger(; kwargs...) # use ProgressMeter for Jupyter    
    elseif isa(stderr, Base.TTY) && (get(ENV, "CI", nothing) ≠ true) # running in terminal?
        global_logger(TerminalLogger()) # enable progress bar
    else
        global_logger(ConsoleLogger())
    end
    Base.global_logger()
end

end
