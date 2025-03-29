dt = 0.01
τ =3
D = 0.01
σY = 0.05
σx = sqrt(D*τ)
α = 0.5
θ_RG = σx,α,σY,dt,0.5
init = [0,1.,1]


function φ(y, y0, σY, α)
    μ = log(2) + α * y0
    A = exp(-((y - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
    B = (1 - erf((y - μ) / (sqrt(2) * σY)))
    return A / B
end

data_all = []
num_lins = 1
num_cells = 50000
arange = collect(-1:0.1:3)
for i in 1:length(arange)
    a = arange[i]
    β(gr,y,y0) = gr^a*φ(y,y0,σY,α)
    θ_OU = 1,D,β,dt
    data= make_lineages(generator_RG,init,θ_RG,num_cells,num_lins);
    df = DataFrame(data,["time","x","length","cell", "lineage"]);
    df = df[df.cell .> 10,:];
    df.time = df.time .- df.time[1];
    df.cell = df.cell .- minimum(df.cell);
    df[:,:T] = cumsum(df.x)*mean(diff(df.time));
    df[:,:lineage] = ones(length(df.time))*i 
    df[:,:a] = a*ones(length(df.time))
    push!(data_all,df)
end
data_OU_scaled = vcat(data_all...);
CSV.write("./data/fig3_data.csv",data_OU_scaled)