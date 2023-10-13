using PhylogeneticTrees
using Test

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
        @test distances[3,5] == 1
        @test distances[6,7] == 1
        @test distances[2,3] == 1
        @test distances[2,5] == 2
        @test distances[3,6] == 3
        @test (1,2) ∉ keys(distances)
        nodes = [2,3,4,5,6,7]
        for i in 1:length(nodes)-1
            for j in i+1:length(nodes)
                @test distances[nodes[i],nodes[j]] == distances[nodes[j],nodes[i]]
            end
        end
    end

    @testset "Disconnected" begin
        tree = PhylogeneticTree([1, 2])
        add_child!(tree, 1, 3)
        add_child!(tree, 2, 4)
        mrca, distances, mrca_distances = compute_pairwise_distances!(tree, Set([3,4]))
        @test (3,4) ∉ keys(distances)
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
            compute_pairwise_distances!(tree, Set([12, 16, 19]), remove_unreachable_nodes=true)
        for key in [1, 2, 4, 5, 6, 10, 11, 12, 16, 19]
            @test key ∈ keys(tree.tree)
        end
        for key in [3, 7, 8, 9, 13, 14, 15, 17, 18]
            @test key ∉ keys(tree.tree)
        end

    end

end
