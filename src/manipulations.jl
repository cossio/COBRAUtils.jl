import COBRA

export reactionsIrrev, addFluxCost


"""
Creates an extended model, where reactions are split in two, 
forward and backward with non-negative fluxes. For each
reaction in the original model, a new reaction is appended
at the end which represents the backward flux.
"""
function reactionsIrrev(cobra_model::COBRA.LPproblem)
    S = [(cobra_model.S) (-cobra_model.S)]
    lb = [max.(0, cobra_model.lb); max.(0, -cobra_model.ub)]
    ub = [max.(0, cobra_model.ub); max.(0, -cobra_model.lb)]
    b = copy(cobra_model.b)
    c = [cobra_model.c; -cobra_model.c]
    rxns = [cobra_model.rxns; cobra_model.rxns .* "_bk"]
    COBRA.LPproblem(S, b, c, lb, ub, cobra_model.osense,
                    copy(cobra_model.csense), rxns,
                    copy(cobra_model.mets))
end


"""
Constrains the sum of all fluxes (weighted by costs) to be less than 1.
Adds a pseudometabolite that is produced by all reactions in the original
model, and a pseudoreaction that is the only consumer of this pseudometabolite.
Then limits the flux of the pseudoreaction to be less than 1.
"""
function addFluxCost(cobra_model::COBRA.LPproblem, costs::Vector{Float64})
    (m,n) = size(cobra_model.S)
    @assert length(costs) == n
    S = [(cobra_model.S) zeros(m); costs' -1]
    lb = [cobra_model.lb; 0]
    ub = [cobra_model.ub; 1]
    b = [cobra_model.b; 0]
    c = [cobra_model.c; 0]
    rxns = [cobra_model.rxns; "PseudoRxnCost"]
    mets = [cobra_model.mets; "PseudoMetCost"]
    csense = [cobra_model.csense; '=']
    COBRA.LPproblem(S, b, c, lb, ub, 
                    cobra_model.osense, csense, 
                    rxns, mets)
end
