using DataFrames, PythonPlot, Distributions,CSV,Random
using SpecialFunctions 
PythonPlot.svg(true)
include("simulations.jl")
include("branching.jl")
include("drawing.jl")
include("single_cell_models.jl")


# FIG 2 ------------------------------------------------------------------------
dt = 0.01
τ =3
D = 0.01
σY = 0.05
σx = sqrt(D*τ)
α = 0.5
θ_RG = σx,α,σY,dt,0.5
init = [0,1.,1]


#include("fig2_generate_lineage_data.jl")
#include("fig2_generate_population_data.jl")

# FIG 3 ------------------------------------------------------------------------
include("fig3_generate_lineage_data.jl")