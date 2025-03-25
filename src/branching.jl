"""
    Cell{T}

A mutable struct representing a cell in a tree.

## Fields

- `label::T`: Data associated with the cell.
- `left::Union{Cell{T}, Nothing}`: Left child cell.
- `right::Union{Cell{T}, Nothing}`: Right child cell.

"""
mutable struct Cell{T}
    label::T
    left::Union{Cell{T}, Nothing}
    right::Union{Cell{T}, Nothing}
end



"""
    create_cell(label, left=nothing, right=nothing)

Create a new cell with the given label and optional left and right children.

## Arguments

- `label`: Data to be stored in the cell.
- `left=nothing`: Left child cell (default is `nothing`).
- `right=nothing`: Right child cell (default is `nothing`).

## Returns

- `Cell`: A new `Cell` object with the specified label and children.

"""
function create_cell(label, left=nothing, right=nothing)
    return Cell(label, left, right)
end

"""
    grow_tree!(cell, terminate, params, generator)

Recursively grow a tree from a single root cell.

## Arguments

- `cell`: The root cell from which to grow the tree.
- `terminate`: A function to determine when to terminate growth.
- `params`: Parameters used in the growth process.
- `generator`: A function to generate new cell labels.

"""
function grow_tree!(cell, terminate, params, generator)
    L = generator(cell.label[end, :], params)
    R = generator(cell.label[end, :], params)
    cell.left = create_cell(L)
    cell.right = create_cell(R)
    if terminate(cell) == false
        grow_tree!(cell.left, terminate, params, generator)
        grow_tree!(cell.right, terminate, params, generator)
    end
end

"""
    grow_forest!(nodes, terminate, params, dist)

Grow a forest of trees from a collection of root nodes.

## Arguments

- `nodes`: Collection of root nodes for trees.
- `terminate`: A function to determine when to terminate growth.
- `params`: Parameters used in the growth process.
- `dist`: Distribution used for generating new cell labels.

"""
function grow_forest!(nodes, terminate, params, generator)
    for node in nodes
        grow_tree!(node, terminate, params, generator)
    end
end

function  make_forest_data_frame(roots::Vector{Cell{T}}) where T
    dfs = [make_tree_data_frame(cell) for cell in roots]
    for (i, df) in enumerate(dfs)
        df.source = fill(i, nrow(df))
    end
    return vcat(dfs...)
end


function sample_or_all(lst, N)
    M = length(lst)
    return N ≤ M ? lst[randperm(M)[1:N]] : lst
end



function grow_dilute!(initial_cells,Td,Tmax,Nmax,params,generator)
    t = 0.0
    current_population = initial_cells


    while t < Tmax
        function terminate(cell)
            cell.label[end,1]>t + Td
        end
        for cell in current_population
            grow_tree!(cell, terminate, params, generator)
        end
        current_population = vcat([get_leaf_nodes(cell) for cell in current_population]...)
        current_population = sample_or_all(current_population, Nmax)
        t += Td
    end
end
"""
    find_max(node::Cell)

Recursively find the maximum value in a tree of cells.

## Arguments

- `node::Cell`: The root node of the tree.

## Returns

- `T`: The maximum value in the tree.

"""
function find_max(node::Cell)
    if node.left != nothing
        left_max = find_max(node.left)
        right_max = find_max(node.right)
        current_max = max(node.label, max(left_max, right_max))
        return current_max
    else
        return node.label
    end
end

"""
    get_leaf_nodes(node)

Get a list of leaf nodes in a tree.

## Arguments

- `node`: The root node of the tree.

## Returns

- `Array`: An array of leaf nodes

"""
function get_leaf_nodes(node)
    if (node.left != nothing) && (node.right != nothing)
        L = get_leaf_nodes(node.left)
        R = get_leaf_nodes(node.right)
        return vcat(L, R)
    else
        return [node]
    end
end

"""
    get_lineage(node)

Get the lineage of a node in the tree.

## Arguments

- `node`: The node for which to retrieve the lineage.

## Returns

- `Array`: An array of node labels representing the lineage.

"""
function get_lineage(node)
    lineage = [node.label]
    while (node.left != nothing) && (node.right != nothing)
        r = rand()
        if r < 0.5
            node = node.left
        else
            node = node.right
        end
        push!(lineage, node.label)
    end
    return vcat(lineage...)
end

"""
    get_node(root, path)

Get a node from the tree given a path.

## Arguments

- `root`: The root node of the tree.
- `path`: An array of 0s and 1s representing the path to the desired node.

## Returns

- `T`: The label of the node at the specified path.

"""
function get_node(root, path)
    node = root
    for p in path[0:-1]
        if p == 0
            node = node.left
        elseif p == 1
            node = node.right
        end
    end
    return node.label
end

function make_forest_data_frame(roots::Vector{Cell{T}}) where T
    dfs = [make_tree_data_frame(cell) for cell in roots]
    for (i, df) in enumerate(dfs)
        df[!, :source] = fill(i, nrow(df))
    end
    return vcat(dfs...)
end

function _collect_tree_data!(node::Cell, branch::String, rows, branch_map::Dict{String, Int}, source::Int)
    branch_id = get!(branch_map, branch) do
        length(branch_map) + 1
    end
    isleaf = (node.left === nothing) && (node.right === nothing)
    for i in 1:size(node.label, 1)
        time, x, length = node.label[i, :]
        push!(rows, (time=time, length=length, x=x, branch=branch, branch_id=branch_id, isleaf=isleaf, source=source))
    end
    if node.left !== nothing
        _collect_tree_data!(node.left, branch * "0", rows, branch_map, source)
    end
    if node.right !== nothing
        _collect_tree_data!(node.right, branch * "1", rows, branch_map, source)
    end
end


function make_tree_data_frame(root::Cell)
    rows = Vector{NamedTuple{(:time, :length, :x, :branch, :branch_id, :isleaf), 
        Tuple{Float64, Float64, Float64, String, Int, Bool}}}()
    branch_map = Dict{String, Int}()
    _collect_tree_data!(root, "0", rows, branch_map)
    return DataFrame(rows)
end

function _collect_tree_data!(node::Cell, branch::String, rows, branch_map::Dict{String, Int})
    branch_id = get!(branch_map, branch) do
        length(branch_map) + 1
    end
    isleaf = (node.left === nothing) && (node.right === nothing)
    for i in 1:size(node.label, 1)
        time, x, length = node.label[i, :]
        push!(rows, (time=time, length=length, x=x, branch=branch, branch_id=branch_id, isleaf=isleaf))
    end
    if node.left !== nothing
        _collect_tree_data!(node.left, branch * "0", rows, branch_map)
    end
    if node.right !== nothing
        _collect_tree_data!(node.right, branch * "1", rows, branch_map)
    end
end
