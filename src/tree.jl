using DataStructures
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
    @assert child_id > parent_id "Child node must have a larger ID than parent"
    @assert parent_id ∈ keys(tree.tree) "Parent node must be in tree"
    @assert child_id ∉ keys(tree.tree) "Child node must not be in tree"
    parent = tree.tree[parent_id]
    child = PhylogeneticNode(child_id, parent, [])
    push!(parent.children, child)   # add child to parent's children
    tree.tree[child_id] = child     # add child to tree
    delete!(tree.leaves, parent_id) # remove parent from leaves
    tree.leaves[child_id] = child   # add child to leaves
end

function compute_pairwise_distances!(tree::PhylogeneticTree, ids::Set{Int}; remove_unreachable_nodes=false)
    """Compute pairwise distances between nodes in the tree, and between 
    nodes in the tree and the MRCA. Distances are computed between nodes in
    `ids` and their ancestors, all the way up to the MRCA.
    

    Params:
        tree::PhylogeneticTree: the tree to compute distances for
        ids::Set{Int}: the set of IDs to start working up the tree

    Returns:
        mrca::Union{Int, Nothing}: the MRCA of the tree
        pairwise_distances::Dict{Tuple{Int,Int}, Int} where pairwise_distances[id1, id2]
            is the distance between id1 and id2
        mrca_distances::Dict{Int, Int} where mrca_distances[id] is the distance between
            id and the MRCA
    """

    # offspring_distances[id1] is the distance between id1 and all of its offspring
    # starting from the ids, we work our way up the tree, merging and 
    # incrementing offspring_distances as we go
    offspring_distances = Dict{Int, Dict{Int, Int}}()

    # pairwise_distances[id1, id2] is the distance between id1 and id2
    # for each (id1, id2) pair, we also store (id2, id1)
    pairwise_distances = Dict{Tuple{Int, Int}, Int}()

    # We use a priority queue to process nodes from oldest to youngest. This
    # guarantees that we will have computed the pairwise distances within the
    # subtree rooted at a node before we process that node.
    pq = PriorityQueue{PhylogeneticNode, Int}(Base.Order.Reverse)

    mrca = nothing
    genesis_nodes = Set{PhylogeneticNode}()

    # start off with leaves/members of pop
    for id in ids
        pq[tree.tree[id]] = id
    end
    # process nodes from oldest to youngest, computing pairwise distances between all
    # children each time as we work our way up the tree
    while length(pq) > 0
        node = dequeue!(pq)

        # node is the MRCA if there are no more nodes in the queue
        # and we haven't seen any other genesis nodes
        if length(pq) == 0 && length(genesis_nodes) == 0
            mrca = node.id
        end
        # add parent if it exists and node is not the MRCA
        if !isnothing(node.parent) && length(pq) > 0
            pq[node.parent] = node.parent.id
        end
        if isnothing(node.parent)
            push!(genesis_nodes, node)
        end

        # compute distance between each offspring and node
        offspring_distances[node.id] = Dict{Int, Int}(node.id => 0)
        pairwise_distances[node.id, node.id] = 0
        for child in node.children
            child.id ∉ keys(offspring_distances) && continue
            for (o_id, o_dist) in offspring_distances[child.id]
                offspring_distances[node.id][o_id] = o_dist + 1
                pairwise_distances[node.id, o_id] = o_dist + 1
                pairwise_distances[o_id, node.id] = o_dist + 1

            end
        end
        # compute pairwise distances between all offspring
        # first two loops: iterate over all pairs of children
        for (cidx1, child1) in enumerate(node.children[1:end-1])
            child1.id ∉ keys(offspring_distances) && continue
            for (cidx2, child2) in enumerate(node.children[cidx1+1:end])
                child2.id ∉ keys(offspring_distances) && continue
                # second two loops: iterate over all pairs of offspring between the two children
                for (o_id1, o_dist1) in offspring_distances[child1.id]
                    for (o_id2, o_dist2) in offspring_distances[child2.id]
                        pairwise_distances[o_id2, o_id1] = o_dist1 + o_dist2 + 2
                        pairwise_distances[o_id1, o_id2] = o_dist1 + o_dist2 + 2
                    end
                end
            end
        end
    end

    # remove unreachable nodes from the tree
    if remove_unreachable_nodes
        for id in keys(tree.tree)
            if id ∉ keys(offspring_distances)
                delete!(tree.tree, id)
                delete!(tree.leaves, id)
            end
        end
    end
    mrca_distances = isnothing(mrca) ? Dict{Int, Int}() : offspring_distances[mrca]
    return mrca, pairwise_distances, mrca_distances
end
