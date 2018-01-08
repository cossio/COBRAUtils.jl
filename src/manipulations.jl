import COBRA

export reactionsIrrev, addFluxCost, remove_reaction, remove_reactions


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


"""
Index of reaction rxn.
"""
function reaction_index(cobra_model::COBRA.LPproblem, rxn::String)
    k = findfirst(cobra_model.rxns, rxn)
    k == 0 && throw(KeyError(rxn))
    return k
end


"""
Index of metabolite met.
"""
function metabolite_index(cobra_model::COBRA.LPproblem, met::String)
    i = findfirst(cobra_model.mets, met)
    i == 0 && throw(KeyError(met))
    return i
end


"""
Returns list of metabolite indices that participate in the k'th reaction.
"""
function reaction_metabolites(cobra_model::COBRA.LPproblem, k::Int)
    @assert 1 ≤ k ≤ length(cobra_model.rxns)
    find(cobra_model.S[1:end, k])
end


"""
Returns list of metabolite indices that are consumed in the k'th reaction.
"""
function reaction_substrates(cobra_model::COBRA.LPproblem, k::Int)
    @assert 1 ≤ k ≤ length(cobra_model.rxns)
    find(s -> s < 0, cobra_model.S[1:end, k])
end


"""
Returns list of metabolite indices that are produced in the k'th reaction.
"""
function reaction_products(cobra_model::COBRA.LPproblem, k::Int)
    @assert 1 ≤ k ≤ length(cobra_model.rxns)
    find(s -> s > 0, cobra_model.S[1:end, k])
end


reaction_metabolites(cobra_model::COBRA.LPproblem, rxn::String) = reaction_metabolites(cobra_model, reaction_index(cobra_model, rxn))
reaction_substrates(cobra_model::COBRA.LPproblem, rxn::String)  = reaction_substrates(cobra_model,  reaction_index(cobra_model, rxn))
reaction_products(cobra_model::COBRA.LPproblem, rxn::String)    = reaction_products(cobra_model,    reaction_index(cobra_model, rxn))


function metabolite_reactions(cobra_model::COBRA.LPproblem, i::Int)
    @assert 1 ≤ i ≤ length(cobra_model.mets)
    find(cobra_model.S[i,:])
end

metabolite_reactions(cobra_model::COBRA.LPproblem, met::String) = metabolite_reactions(cobra_model, metabolite_index(cobra_model, met))


"""
Returns list of reactions that exchange the i'th metabolite
"""
function exchange_reactions(cobra_model::COBRA.LPproblem, i::Int)
    indexes = Int[]
    for k in metabolite_reactions(cobra_model, i)
        if isempty(reaction_substrates(cobra_model, k)) || isempty(reaction_products(cobra_model, k))
            push!(indexes, k)
        end
    end
    return indexes
end

exchange_reactions(cobra_model::COBRA.LPproblem, met::String) = exchange_reactions(cobra_model, metabolite_index(cobra_model, met))


"""
Returns a copy of cobra_model where reaction i has been eliminated
"""
remove_reaction(cobra_model::COBRA.LPproblem, i::Int) = remove_reactions(cobra_model, [i])


"""
Returns a copy of cobra_model where the reactions in list have been eliminated
"""
function remove_reactions(cobra_model::COBRA.LPproblem, list::AbstractVector{Int})
    keep = [k for k = 1 : length(cobra_model.rxns) if k ∉ list]
    keep_reactions(cobra_model, keep)
end


"""
Returns a copy of cobra_model where all reactions have been eliminated except those in list.
"""
function keep_reactions(cobra_model::COBRA.LPproblem, list::AbstractVector{Int})
    @assert all(1 .≤ list .≤ length(cobra_model.rxns))
    list = unique(list)

    S = cobra_model.S[:, list]
    c = cobra_model.c[list]
    lb = cobra_model.lb[list]
    ub = cobra_model.ub[list]
    rxns = cobra_model.rxns[list]

    COBRA.LPproblem(S, cobra_model.b, c, lb, ub, cobra_model.osense, 
                    cobra_model.csense, rxns, cobra_model.mets)
end


"""
Returns a copy of model where metabolites in list have been eliminated.
"""
function remove_metabolites(model::COBRA.LPproblem, list::AbstractVector{Int})
    keep = [i for i = 1 : length(model.mets) if i ∉ list]
    keep_metabolites(model, keep)
end


"""
Returns a copy of model where metabolites in list have been eliminated.
"""
remove_metabolite(model::COBRA.LPproblem, i::Int) = remove_metabolites(model, [i])


"""
Returns a copy of model where all metabolites have been eliminated except those in list.
"""
function keep_metabolites(model::COBRA.LPproblem, list::AbstractVector{Int})
    @assert all(1 .≤ list .≤ length(model.mets))
    list = unique(list)

    S = model.S[list, :]
    b = model.b[list]
    csense = model.csense[list]
    mets = model.mets[list]

    COBRA.LPproblem(S, b, model.c, model.lb, model.ub, 
                    model.osense, csense, model.rxns, mets)
end


"""
Returns a copy of model where reactions that can't carry flux have been removed.
Null reactions are those for which max(|fmin|, |fmax|) < tol, where tol = 1e-4
by default. This requires fmin, fmax from fva. If you have computed fmin, fmax 
already, you can pass them as arguments to save computing time.
"""
function reduce_model end

function reduce_model(model::COBRA.LPproblem; tol::Real = 1e-6)
    @assert 0 ≤ tol < Inf
    fmin, fmax = fva(model)
    reduce_model(model, fmin, fmax; tol = tol)
end

function reduce_model(model::COBRA.LPproblem, fmin::Vector{Float64}, fmax::Vector{Float64}; 
                      tol::Real = 1e-6)
    @assert 0 ≤ tol < Inf
    @assert length(fmin) == length(fmax) == length(model.rxns)
    rxns = [k for k = 1 : length(model.rxns) if abs(fmin[k]) < tol && abs(fmax[k]) < tol]
    model = remove_reactions(model, rxns)
    mets = [i for i = 1 : length(model.mets) if length(metabolite_reactions(model, i)) == 0]
    model = remove_metabolites(model, mets)
    return model
end


"""
Sets the objective reaction as the i'th reaction.
If sense == 1 (default), flux is maximized. If
sense == -1, flux is minimized.
"""
function set_objective!(model::COBRA.LPproblem, i::Int, sense::Int = 1)
    model.c .= 0
    model.c[i] = 1
    model.osense = sign(sense)
end