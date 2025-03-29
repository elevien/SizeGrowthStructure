# setup model parameters --------------------------------------------------------
dt = 0.01
τ = 1
D = 0.02
σY = 0.05
σx = sqrt(D*τ)
α = 0.5
init = [0,1.,1]

function φ(y, y0, σY, α)
    μ = log(2) + α * y0
    A = exp(-((y - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
    B = (1 - erf((y - μ) / (sqrt(2) * σY)))
    return A / B
end

# Define the base parameters
θ_OU_rate = (τ, D, (gr,y,y0) -> gr*φ(y,y0,σY,α), dt)

# generate seed population ------------------------------------------------------
Tmax_init = 5
function terminate(cell)
    cell.label[end,1]>Tmax_init
end

root = create_cell(generator_OU_rate([-Tmax_init,1.,1],θ_OU_rate))
grow_tree!(root,terminate,θ_OU_rate,generator_OU_rate)
initial_population = get_leaf_nodes(root)

# loop over arange --------------------------------------------------------------
Tmax = 200
Td = 1
Nmax = 800

data_all = []
arange = collect(-1:0.1:3)
for i in 1:length(arange)
    # Create a fresh copy of the initial population for each a value
    population = deepcopy(initial_population)
    a = arange[i]
    println("running a = $a")
    
    # Create parameters with current a value
    β(gr,y,y0) = gr^a*φ(y,y0,σY,α)
    θ_OU = (τ, D, β, dt)
    
    grow_dilute!(population, Td, Tmax, Nmax, θ_OU, generator_OU_rate)
    df = make_forest_data_frame(population)
    df[!, :time_r] = round.(df.time, digits=1)
    df = df[df.time_r .< Tmax,:]
    df[!, :a] = fill(a, nrow(df))
    df = df[df.time_r .>= maximum(df.time_r),:]
    push!(data_all,df)
end
data_all = vcat(data_all...)
CSV.write("./data/fig2_population_data.csv",data_all)