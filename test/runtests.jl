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
        mrca, distances, mrca_distances = compute_pairwise_distances(tree, Set([5,7]))
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
        mrca, distances, mrca_distances = compute_pairwise_distances(tree, Set([3,4]))
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
        mrca, distances, mrca_distances = compute_pairwise_distances(tree, Set([3,4]))
        @test mrca == 1
        add_child!(tree, 3, 5)
        add_child!(tree, 4, 6)
        mrca, distances, mrca_distances = compute_pairwise_distances(tree, Set([5,6]))
        @test mrca == 1
    end

    @testset "Compute distance between subset of leaves" begin
        # This test makes sure that we don't need to process 4 before processing 1
        tree = PhylogeneticTree([1, 2])
        add_child!(tree, 1, 3)
        add_child!(tree, 1, 4)
        mrca, distances, mrca_distances = compute_pairwise_distances(tree, Set([2, 3]))
        @test isnothing(mrca)
        @test (1,2) ∉ keys(distances)
        @test (1,4) ∉ keys(distances)
        @test distances[1,3] == 1
    end

    @testset "Large Random Tree, all pairwise distances" begin
        tree = PhylogeneticTree(collect(1:100))
        for i in 101:1000
            add_child!(tree, i - rand(1:10), i+1)
        end
        ids = Set(collect(900:1000)) 
        mrca, distances, mrca_distances = compute_pairwise_distances(tree, ids)
        for id1 in ids
            for id2 in ids
                id1 == id2 && continue
                @test distances[id1, id2] == distances[id2, id1]
            end
        end
    end


end
