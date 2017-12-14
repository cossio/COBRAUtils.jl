import COBRA

export AVCobraModel


"""
Minimal metabolic network exhibiting a switch, in COBRA's format.
The network objective is ATP production.
    Vg, Vl: uptake bounds of glucose and lactate
    R: maximum respiratory rate
    L: large bound, used in place of ∞ (default: 1000)
"""
function AVCobraModel(Vg::Real, Vl::Real = 0; R = 0.9, L = 1000)
    @assert 0 ≤ Vg < Inf && 0 ≤ Vl < Inf
    @assert 0 ≤ R < Inf
    @assert 0 < L < Inf

    rxns = ["vg", "vl", "vo", "vatp"]
    mets = ["pyr", "atp"]

    S = [2 1 -1 0; 2 0 18 -1]
    b = zeros(2)
    c = [0, 0, 0, 1]
    lb = [-L, -L, 0, -L]
    ub = [Vg, Vl, R, L]
    osense = 1
    csense = ['E', 'E']

    COBRA.LPproblem(S, b, c, lb, ub, osense, csense, rxns, mets)
end
