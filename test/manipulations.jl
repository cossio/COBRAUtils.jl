using Base.Test
import COBRAUtils

@testset "manipulations" begin


@testset "remove reaction" begin
    model = COBRAUtils.loadEColiTestModel()
    n = length(model.rxns)
    model2 = COBRAUtils.remove_reaction(model, 1)
    @test length(model2.rxns) == length(model.rxns) - 1
    @test length(model2.lb) == length(model.lb) - 1
    @test length(model2.ub) == length(model.ub) - 1
    @test length(model2.c) == length(model.c) - 1
    @test size(model2.S) == (size(model.S, 1), n - 1)
end


@testset "remove reaction list" begin
    model = COBRAUtils.loadEColiTestModel()
    n = length(model.rxns)
    model2 = COBRAUtils.remove_reactions(model, [1,2])
    @test length(model2.rxns) == length(model.rxns) - 2
    @test length(model2.lb) == length(model.lb) - 2
    @test length(model2.ub) == length(model.ub) - 2
    @test length(model2.c) == length(model.c) - 2
    @test size(model2.S) == (size(model.S, 1), n - 2)
end


end

