import COBRA


function lpproblem(S::AbstractMatrix, b::AbstractVector, c::AbstractVector, 
                   lb::AbstractVector, ub::AbstractVector)
    @assert length(b) == size(S,1)
    @assert length(lb) == length(ub) == length(c) == size(S,2)
    osense = -1 # maximize
    csense = fill('=', length(b))
    mets = "M" .* string.(1 : size(S,1))
    rxns = "R" .* string.(1 : size(S,2))
    COBRA.LPproblem(S, b, c, lb, ub, osense, csense, rxns, mets)
end


"""
Converts a MetNet to COBRA.LPProblem.
The MetNet struct is defined in
https://github.com/anna-pa-m/Metabolic-EP
"""
function lpproblem(metnet)
    lpproblem(metnet.S, metnet.b, metnet.c, metnet.lb, metnet.ub)
end