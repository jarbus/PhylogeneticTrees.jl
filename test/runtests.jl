using PhylogeneticTrees
using Test
using Serialization

@testset "Construction" begin
    n_pop = 10
    genesis_pop = collect(1:n_pop)
    tree = PhylogeneticTree(genesis_pop)
    @testset "Genesis" begin
        @test length(tree.tree) == n_pop
        @test length(tree.genesis) == n_pop
        for gen_ind in tree.genesis
            @test tree.tree[gen_ind.id] == gen_ind
        end
        @test isnothing(tree.mrca)
        @test tree.tree == tree.leaves
    end
    @testset "add_child!" begin
        parent_ids = collect(1:n_pop)
        child_ids = collect(n_pop+1:2n_pop)
        for (parent_id, child_id) in zip(parent_ids, child_ids)
            add_child!(tree, parent_id, child_id)
        end

        for (parent_id, child_id) in zip(parent_ids, child_ids)
            @test tree.tree[child_id].id == child_id
            @test tree.tree[child_id].parent.id == parent_id
            @test child_id in [c.id for c in tree.tree[child_id].parent.children]
            @test tree.leaves[child_id] == tree.tree[child_id]
            @test tree.leaves[child_id].parent.id ∉ keys(tree.leaves)
        end

        @test length(tree.tree) == 2n_pop
        @test length(tree.genesis) == n_pop
        @test length(tree.leaves) == n_pop
    end
end

@testset "Distances" begin
    @testset "Connected" begin
        #   1
        #   |
        #   2
        #  / \
        # 3   4
        # |   |
        # 5   6
        #     |
        #     7
        tree = PhylogeneticTree([1])
        add_child!(tree, 1, 2)
        add_child!(tree, 2, 3)
        add_child!(tree, 2, 4)
        add_child!(tree, 3, 5)
        add_child!(tree, 4, 6)
        add_child!(tree, 6, 7)
        mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set([5,7]))
        @test distances[5,7] == 5
        @test distances[3,4] == 2
        @test distances[3,7] == 4
        @test distances[3,3] == 0
        @test distances[3,5] == 1
        @test distances[6,7] == 1
        @test distances[2,3] == 1
        @test distances[2,5] == 2
        @test distances[3,6] == 3
        @test (1,2) ∉ keys(distances)
        nodes = [2,3,4,5,6,7]
    end

    @testset "Disconnected" begin
        tree = PhylogeneticTree([1, 2])
        add_child!(tree, 1, 3)
        add_child!(tree, 2, 4)
        mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set([3,4]))
        @test (3,4) ∉ keys(distances)
        @test distances[1,1] == 0
        @test distances[1,3] == 1
        @test distances[2,4] == 1
        @test isnothing(mrca)
        @test length(mrca_distances) == 0
    end

    @testset "Connected MRCA at genesis" begin
        tree = PhylogeneticTree([1, 2])
        add_child!(tree, 1, 3)
        add_child!(tree, 1, 4)
        mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set([3,4]))
        @test mrca == 1
        add_child!(tree, 3, 5)
        add_child!(tree, 4, 6)
        mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set([5,6]))
        @test mrca == 1
    end

    @testset "Compute distance between subset of leaves" begin
        # This test makes sure that we don't need to process 4 before processing 1
        # Also checks that we can remove unreachable nodes if specified
        #    1   2
        #   / \
        #  3   4
        tree = PhylogeneticTree([1, 2])
        add_child!(tree, 1, 3)
        add_child!(tree, 1, 4)
        mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set([2, 3]))
        @test isnothing(mrca)
        @test (1,2) ∉ keys(distances)
        @test (1,4) ∉ keys(distances)
        @test distances[1,3] == 1
        mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set([2, 3]), remove_unreachable_nodes=true)
        @test 4 ∉ keys(tree.tree) && 4 ∉ keys(tree.leaves)
        @test 1 ∈ keys(tree.tree)
        @test 2 ∈ keys(tree.tree)
    end

    @testset "Remove unreachable nodes" begin
        #     1          2      3
        #    / \        / \    / \
        #   4   5      6   7  8   9
        #  /     \     |   |  |   |
        # 10     11   12  13 14  15
        # / \    / \
        #16 17  18 19
        tree = PhylogeneticTree([1, 2, 3])
        add_child!(tree, 1, 4)
        add_child!(tree, 1, 5)
        add_child!(tree, 2, 6)
        add_child!(tree, 2, 7)
        add_child!(tree, 3, 8)
        add_child!(tree, 3, 9)
        add_child!(tree, 4, 10)
        add_child!(tree, 5, 11)
        add_child!(tree, 6, 12)
        add_child!(tree, 7, 13)
        add_child!(tree, 8, 14)
        add_child!(tree, 9, 15)
        add_child!(tree, 10, 16)
        add_child!(tree, 10, 17)
        add_child!(tree, 11, 18)
        add_child!(tree, 11, 19)
        mrca, distances, mrca_distances =
            compute_pairwise_distances!(tree, Set([12, 16, 19]),
                                        remove_unreachable_nodes=true,
                                       max_distance=2)
        for key in [1, 2, 4, 5, 6, 10, 11, 12, 16, 19]
            @test key ∈ keys(tree.tree)
        end
        for key in [3, 7, 8, 9, 13, 14, 15, 17, 18]
            @test key ∉ keys(tree.tree)
        end

    end

end

@testset "Perf" begin
    # create a tree with 8192 leaves, each node has two children
    function recursively_add_two_children(tree, id, depth)
        if depth == 0
            return
        end
        add_child!(tree, id, id*2+1)
        add_child!(tree, id, id*2+2)
        recursively_add_two_children(tree, id*2+1, depth-1)
        recursively_add_two_children(tree, id*2+2, depth-1)
    end

    tree = PhylogeneticTree([0])
    depth = 13
    recursively_add_two_children(tree, 0, depth)
    max_distance = 12
    mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set(collect(2^(depth-1)+1:2^depth)), max_distance=max_distance)
end

@testset "Memory" begin
    tree = PhylogeneticTree(collect(1:100))
    n_gens = 1000
    max_distance = 5
    for gen in 1:n_gens
        # Print memory usage
        parent_ids = rand(  (100(gen-1))+1:(100gen), 10)
        child_ids = collect((100gen)+1  :((100gen)+100))
        # Add children randomly across the 10 parents
        for child_id in child_ids
            add_child!(tree, rand(parent_ids), child_id)
        end
        compute_pairwise_distances!(tree, 
                                    Set(child_ids),
                                    max_distance=max_distance,
                                    remove_unreachable_nodes=true)
    end
end

@testset "Equality" begin
    # Test that two nodes are equal if they have the same ID, no parents and no children
    node1 = PhylogeneticNode(1, nothing, [])
    node2 = PhylogeneticNode(1, nothing, [])
    node3 = PhylogeneticNode(2, nothing, [])
    node4 = PhylogeneticNode(1, node1, [])
    node5 = PhylogeneticNode(1, nothing, [node1])
    node6 = PhylogeneticNode(1, nothing, [node1, node2])

    node1.id == node2.id || println("a")
    # check that parents are either both nothing or their ids are equal
    @test node1 == node2
    @test node1 != node3
    @test node1 != node4
    @test node5 != node6
    @test node6 == node6
end

@testset "Serialization" begin
   struct TreeContainer
       tree::PhylogeneticTree
   end
   tree_original = PhylogeneticTree([1])
   # add a chain of children
   for i in 2:5_001
       add_child!(tree_original, i-1, i)
   end
   container = TreeContainer(tree_original)
   # Check node serialization
   isfile("t.jls") && rm("t.jls")
   serialize("t.jls", container)
   tree_deserialized = deserialize("t.jls").tree
   # We can test for equality
   @test tree_original.tree[1] == tree_deserialized.tree[1]
   # Add a lot more children
   for i in 5_002:100_000
       add_child!(tree_original, i-1, i)
   end
   # Test serialization and deserialization
   isfile("t.jls") && rm("t.jls")
   serialize("t.jls", container)
   deserialized_container = deserialize("t.jls")
   @test true # we don't test for equality because the tree is too large
end
