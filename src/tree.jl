using DataStructures
using Serialization
export serialize, deserialize
# This needs to be mutable to we can update the references 
# to other nodes for the garbage collector
mutable struct PhylogeneticNode
    id::Int
    parent::Union{PhylogeneticNode, Nothing}
    children::Vector{PhylogeneticNode}
end

function Base.:(==)(node1::PhylogeneticNode, node2::PhylogeneticNode)::Bool
    # Checks that the structure and IDs of all nodes are equal, does not care
    # about the memory location of the nodes
    node1.id == node2.id || return false
    # check that parents are either both nothing or their ids are equal
    isnothing(node1.parent) == isnothing(node2.parent) || return false
    isnothing(node1.parent) || node1.parent.id == node2.parent.id || return false
    # Check that all children are equal
    length(node1.children) == length(node2.children) || return false
    all(i -> node1.children[i] == node2.children[i], 1:length(node1.children))
end


mutable struct PhylogeneticTree
    genesis::Vector{PhylogeneticNode}
    tree::Dict{Int,PhylogeneticNode}
    leaves::Dict{Int,PhylogeneticNode}
    mrca::Union{PhylogeneticNode, Nothing}
end

# Serialization of large trees causes a stack overflow error,
# so we convert the tree to a dict
function Serialization.serialize(s::AbstractSerializer, tree::PhylogeneticTree)
    serialized_tree = Dict{Int, Any}()
    for node in values(tree.tree)
        serialized_tree[node.id] = [isnothing(node.parent) ? 0 : node.parent.id,
                                    [child.id for child in node.children]]
    end
    # We specify that we are writing a PhylogeneticTree so the correct deserialization
    # function is called
    Serialization.writetag(s.io, Serialization.OBJECT_TAG)
    Serialization.serialize(s, PhylogeneticTree)
    Serialization.serialize(s, serialized_tree)
end

function Serialization.deserialize(s::AbstractSerializer, ::Type{PhylogeneticTree})
    serialized_tree = Serialization.deserialize(s)
    # Create empty tree
    tree = PhylogeneticTree(Int[])
    # Add all nodes to tree
    for (id, (parent_id, child_ids)) in serialized_tree
        tree.tree[id] = PhylogeneticNode(id, nothing, [])
    end
    # Add connections between nodes
    for (id, (parent_id, child_ids)) in serialized_tree
        node = tree.tree[id]
        node.parent = parent_id == 0 ? nothing : tree.tree[parent_id]
        node.children = [tree.tree[child_id] for child_id in child_ids if child_id in keys(tree.tree)]
        if length(node.children) == 0
            tree.leaves[id] = node
        end
        if isnothing(node.parent)
            push!(tree.genesis, node)
        end
    end
    return tree
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
    @assert child_id > parent_id "Child node must have a larger ID than parent, but got $child_id <= $parent_id"
    @assert parent_id ∈ keys(tree.tree) "Parent node $parent_id must be in tree"
    @assert child_id ∉ keys(tree.tree) "Child node must not be in tree"
    parent = tree.tree[parent_id]
    child = PhylogeneticNode(child_id, parent, [])
    push!(parent.children, child)   # add child to parent's children
    tree.tree[child_id] = child     # add child to tree
    delete!(tree.leaves, parent_id) # remove parent from leaves
    tree.leaves[child_id] = child   # add child to leaves
end

function compute_pairwise_distances!(tree::PhylogeneticTree,
                                     ids::Set{Int};
                                     remove_unreachable_nodes=false,
                                     max_distance::Int=typemax(Int))
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

    # pairwise_distances[i] is a tuple [id1, id2, dist] where
    # dist is the distance between id1 and id2, and id1 <= id2
    num_nodes = length(tree.tree)
    pairwise_distances = Vector{NTuple{3,Int}}(undef, Int((num_nodes*(num_nodes+1))/2))
    pair_id = 1

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
        pairwise_distances[pair_id] = (node.id, node.id, 0)
        pair_id += 1
        for child in node.children
            child.id ∉ keys(offspring_distances) && continue
            for (o_id, o_dist) in offspring_distances[child.id]
                distance = o_dist + 1
                distance > max_distance && continue
                offspring_distances[node.id][o_id] = distance
                min_node_id = min(node.id, o_id)
                max_node_id = max(node.id, o_id)
                pairwise_distances[pair_id] = (min_node_id, max_node_id, distance)
                pair_id += 1

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
                    o_dist1 > max_distance && continue
                    for (o_id2, o_dist2) in offspring_distances[child2.id]
                        distance = o_dist1 + o_dist2 + 2
                        distance > max_distance && continue
                        min_node_id = min(o_id1, o_id2)
                        max_node_id = max(o_id1, o_id2)
                        pairwise_distances[pair_id] = (min_node_id, max_node_id, distance)
                        pair_id += 1
                    end
                end
            end
        end
    end

    # Convert pairwise_distances to a dictionary at the end to avoid rehashing
    pairwise_distances_dict = Dict{Tuple{Int, Int}, Int}((n1, n2) => d
                        for (n1, n2, d) in pairwise_distances[1:pair_id-1])

    # remove unreachable nodes from the tree
    if remove_unreachable_nodes
        # set parent of mrca to nothing
        # this allows the garbage collector to remove the entire tree from memory
        if !isnothing(mrca) && !isnothing(tree.tree[mrca].parent)
            tree.tree[mrca].parent = nothing
        end
        # remove all genesis nodes that are not in offspring_distances
        tree.genesis = [node for node in tree.genesis if node.id ∈ keys(offspring_distances)]
        # remove all nodes that are not in offspring_distances
        for id in keys(tree.tree)
            if id ∉ keys(offspring_distances)
                delete!(tree.tree, id)
                delete!(tree.leaves, id)
            end
        end
    end
    mrca_distances = isnothing(mrca) ? Dict{Int, Int}() : offspring_distances[mrca]
    return mrca, pairwise_distances_dict, mrca_distances
end
