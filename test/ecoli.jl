using Base.Test
import COBRAUtils


@testset "E. Coli" begin
    model = COBRAUtils.loadEColiTestModel()

    @test COBRAUtils.metabolite_index(model, "6pgc[c]") == 4
    @test COBRAUtils.reaction_index(model, "ACONTa") == 4
    @test sort(model.mets[COBRAUtils.reaction_metabolites(model, "ACALD")]) == sort(["acald[c]", "accoa[c]", "coa[c]", "h[c]", "nad[c]", "nadh[c]"])
    @test sort(model.mets[COBRAUtils.reaction_substrates(model, "ACALD")]) == sort(["acald[c]", "coa[c]", "nad[c]"])
    @test sort(model.mets[COBRAUtils.reaction_products(model, "ACALD")]) == sort(["accoa[c]", "h[c]", "nadh[c]"])
    @test sort(model.rxns[COBRAUtils.metabolite_reactions(model, "ac[c]")]) == sort(["ACt2r", "ACKr"])

    @test COBRAUtils.reaction_substrates(model, "EX_ac(e)") == [7]
    @test isempty(COBRAUtils.reaction_products(model, "EX_ac(e)"))
    
    @test COBRAUtils.exchange_reactions(model, "ac[e]") == [20]

    status, μ, v = COBRAUtils.optimize(model)
    @test μ ≈ 0.8739215069684305
end


@testset "reduce_model E. Coli" begin
    ecoli = COBRAUtils.loadEColiTestModel();
    ecolired = COBRAUtils.reduce_model(ecoli);
    @test size(ecolired.S) == (68, 87)
    status0, μ0, v0 = COBRAUtils.optimize(ecoli)
    status1, μ1, v1 = COBRAUtils.optimize(ecolired)
    @test μ0 ≈ μ1
end
