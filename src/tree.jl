struct PhylogeneticNode
    id::Int
    parent::Union{PhylogeneticNode, Nothing}
    children::Vector{PhylogeneticNode}
end

mutable struct PhylogeneticTree
    genesis::Vector{PhylogeneticNode}
    tree::Dict{Int,PhylogeneticNode}
    leaves::Dict{Int,PhylogeneticNode}
    mrca::Union{PhylogeneticNode, Nothing}
    # distances::Dict{Int, Dict{Int, Int}}
end

function PhylogeneticTree()
    error("Genesis population must be provided for PhylogeneticTree")
end
function PhylogeneticTree(genesis_pop_ids::Vector{Int})
    genesis = [PhylogeneticNode(id, nothing, []) for id in genesis_pop_ids]
    tree = Dict(node.id => node for node in genesis)
    leaves = Dict(node.id => node for node in genesis)
    return PhylogeneticTree(genesis, tree, leaves, nothing)
end

function add_child!(tree::PhylogeneticTree, parent_id::Int, child_id::Int)
    @assert parent_id ∈ keys(tree.tree) "Parent node must be in tree"
    @assert child_id ∉ keys(tree.tree) "Child node must not be in tree"
    parent = tree.tree[parent_id]
    child = PhylogeneticNode(child_id, parent, [])
    push!(parent.children, child)   # add child to parent's children
    tree.tree[child_id] = child     # add child to tree
    delete!(tree.leaves, parent_id) # remove parent from leaves
    tree.leaves[child_id] = child   # add child to leaves
end


# function add_distance!(tree::PhylogeneticTree, id1::Int, id2::Int, distance::Int)
#     if id1 ∉ keys(tree.distances)
#         tree.distances[id1] = Dict{Int, Int}()
#     end
#     tree.distances[id1][id2] = distance
#     tree.distances[id2][id1] = distance
# end
#
# function compute_n_offspring(tree::PhylogeneticTree)
#     n_offspring = Dict{Int, Int}()
#     for (id, node) in tree.leaves
#         n_offspring[id] = 0
#     end
# end
#
# function compute_distances!(tree::PhylogeneticTree)
#     cur_pointers = Dict{PhylogeneticNode, Vector{Vector{Int}}}()
#     next_pointers = Dict{PhylogeneticNode, Vector{Vector{Int}}}()
#     for node in values(tree.leaves)
#         cur_pointers[node] = [[0]]
#     end
#     while length(cur_pointers) > 0
#     end
# end
#
