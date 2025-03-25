using DataFrames, PythonPlot, Distributions,CSV,Random
using SpecialFunctions 
PythonPlot.svg(true)
include("simulations.jl")
include("branching.jl")
include("drawing.jl")
include("single_cell_models.jl")

include("fig2_generate_population_data.jl")