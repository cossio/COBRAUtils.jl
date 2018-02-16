"Describes a component of growth media"
struct MediaComponent
    "metabolite index in a COBRA model"
    metabolite::Int
    "index of corresponding exchange reaction"
    exchange_reaction::Int
    "Maximum uptake rate"
    V::Float64
    "concentration of the metabolite in the media"
    concentration::Float64

    function MediaComponent(i::Integer, ex::Integer, V::Real, c::Real)
        @assert i > 0 && ex > 0 && V ≥ 0 && c ≥ 0
        new(i, ex, V, c)
    end
end

function MediaComponent(cobra_model::COBRA.LPproblem, metabolite::Integer, V::Real, c::Real)
    rxns = exchange_reaction(cobra_model, metabolite)
    @assert length(rxns) == 1
    MediaComponent(metabolite, first(rxns), V, c)    
end


"A growth media is a collection of media components."
Media = Vector{MediaComponent}


"""
Sets uptake bounds: u ≤ min(V, c/ξ)
"""
function set_bounds_ξ(cobra_model::COBRA.LPproblem, media::Media, ξ::Real)
    @assert 0 ≤ ξ < Inf
    for x in media
        if x.c == 0
            cobra_model.ub[x.exchange_reaction] = 0
        elseif ξ == 0
            cobra_model.ub[x.exchange_reaction] = x.V
        else
            cobra_model.ub[x.exchange_reaction] = min(x.V, x.c / ξ)
        end
    end
end
