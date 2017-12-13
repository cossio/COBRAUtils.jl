using Base.Test
import COBRAUtils


function loadEColi()
    path = "../metabolic_networks/ecoli_core_model.mat"
    model = COBRAUtils.readCobraModel(path)
    model.osense = 1
    return model
end


@testset "E. Coli" begin
    model = loadEColi()

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
    @test μ == 0.8739215069684305
end
