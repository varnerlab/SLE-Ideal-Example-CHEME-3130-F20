# Build my environment -
import Pkg
Pkg.activate()
Pkg.instantiate()

# include external packages -
using Plots

# include my code -
include("./src/Compute.jl")