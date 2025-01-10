module WangLandau

using StaticArrays: SVector, MVector
using CommonSolve: CommonSolve, init, solve, solve!#, step!
using Random: Random, seed!, shuffle!
import TOML

const PACKAGE_VERSION = VersionNumber(TOML.parsefile(pkgdir(@__MODULE__, "Project.toml"))["version"])

@doc """
    WangLandau

This is `WangLandau` $PACKAGE_VERSION.
Read the documentation [online]().
"""
WangLandau

export init, solve!, solve  # from CommonSolve
export WangLandauProblem

include("WangLandauProblem.jl")

end
