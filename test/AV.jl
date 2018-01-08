using Base.Test
import COBRAUtils

@testset "AV, 1" begin
    model = COBRAUtils.AVCobraModel(0.5, 0.)
    status, μ, v = COBRAUtils.optimize(model)
    @test status == :Optimal
    @test μ ≈ 17.2
    @test v ≈ [0.5, -0.1, 0.9, 17.2]

    model = COBRAUtils.reactionsIrrev(model)
    status, μ, v = COBRAUtils.optimize(model)
    @test status == :Optimal
    @test μ ≈ 17.2
    @test v ≈ [0.5, 0., 0.9, 17.2, 0.0, 0.1, 0.0, 0.0]

    model = COBRAUtils.addFluxCost(model, ones(length(model.rxns)))
    status, μ, v = COBRAUtils.optimize(model)
    @test status == :Optimal
    @test sum(v[1 : end - 1]) ≈ v[end] ≈ 1
    @test μ ≈ 0.9268292682926829
end


@testset "AV fva" begin
    model = COBRAUtils.AVCobraModel(0.5, 0.)
    fmin, fmax = COBRAUtils.fva(model)
    @test fmin ≈ [0.0, -1.0, 0.0, 0.0] 
    @test fmax ≈ [0.5, 0.0, 0.9, 17.2]
end


@testset "AV find" begin
    model = COBRAUtils.AVCobraModel(0.5, 0.)
    @test COBRAUtils.metabolite_reactions(model, "atp") == [1,3,4]
    @test COBRAUtils.metabolite_reactions(model, "pyr") == [1,2,3]
end