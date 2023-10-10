using PhylogeneticTrees
using Test

@testset "PhylogeneticTrees.jl" begin
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
