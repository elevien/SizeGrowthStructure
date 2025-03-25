function make_lineages(generator,init,theta,num_cells,num_lins)

    data = []
    for lin in 1:num_lins

        X = generator(init,theta)
        X = hcat(X,zeros(length(X[:,1]))) # add column for cell label
        lineage =[X]
        for k in 2:num_cells
            X = generator(X[end,:],theta)
            X = hcat(X,ones(length(X[:,1]))*k) # add column for cell label
            push!(lineage,X)
        end
        lineage = vcat(lineage...)
        lineage = hcat(lineage,ones(length(lineage[:,1]))*lin)
        push!(data,lineage)
    end
    # columns are 
    # time, growth rate,cell, lineage
    data = vcat(data...);
    return data
end

